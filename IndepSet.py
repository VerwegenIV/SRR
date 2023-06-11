from gurobipy import *
import math
import SmallCG
import networkx as nx
import MakeInstance
import pandas as pd
import time


def find_odd_cycle(vardict):
    dict_keys = list(vardict.keys())
    odd_cycle = []
    odd_cycle_weight = math.inf
    G = nx.Graph()
    for key in dict_keys:
        G.add_nodes_from([key[4:-1] + 'e', key[4:-1] + 'o'])
    for i in range(len(dict_keys) - 1):
        var1 = dict_keys[i]
        for j in range(i + 1, len(dict_keys)):
            var2 = dict_keys[j]
            if len(set(vardict[var1]['edges']).intersection(set(vardict[var2]['edges']))) > 0 or vardict[var1]['r'] == \
                    vardict[var2]['r']:
                G.add_edge(var1[4:-1] + 'e', var2[4:-1] + 'o',
                           weight=max(0.00000001, 1 - vardict[var1]['value'] - vardict[var2]['value']))
                G.add_edge(var1[4:-1] + 'o', var2[4:-1] + 'e',
                           weight=max(0.00000001, 1 - vardict[var1]['value'] - vardict[var2]['value']))
    print(list(G.edges))
    for key in dict_keys:
        try:
            weight, path = nx.single_source_dijkstra(G, key[4:-1] + 'e', key[4:-1] + 'o')
            path = set(node[0:-1] for node in path)
        except nx.exception.NetworkXNoPath:
            continue
        if weight < odd_cycle_weight:
            odd_cycle = path
            odd_cycle_weight = weight
            print('found better cycle', odd_cycle, "with weight", weight)
    return odd_cycle, odd_cycle_weight


def solve_dual(model, n, edges, costdict, posvar):
    totalsum = 1

    while totalsum > 0.001:
        # dual vars for rounds constraints
        dual_r = model.Pi[:n - 1]
        # dual vars for match constraints
        dual_m = model.Pi[n - 1:n + int((n - 1) * n / 2) - 1]

        for r in range(n - 1):
            print('Trying round', r)
            maxdual = Model("dual")
            maxdual.ModelSense = GRB.MAXIMIZE
            maxdual.Params.LogToConsole = 0

            x = maxdual.addVars(range(len(edges)), vtype=GRB.BINARY, name='x')
            v = maxdual.addVars(posvar, vtype=GRB.INTEGER, name='v')
            # Objective of dual_r[r] + sum_m(dual_m[m] - cost[m])*x[m]
            maxdual.setObjective(dual_r[r] + quicksum((dual_m[e] - costdict.get((edges[e][0], edges[e][1], r), 0))*x[e]
                                                      for e in range(len(edges))))
            # Constraint to make sure a perfect matching is chosen
            maxdual.addConstrs((quicksum(x[e] for e in range(len(edges)) if v in edges[e]) == 1 for v in range(n)),
                               'matching')
            maxdual.addConstrs((v[key] - quicksum(x[e] for e in range(len(edges)) if edges[e] in posvar[key]['edges'])
                                == 0 for key in posvar), 'AND')
            maxdual.addConstrs((v[key] <= n/2 - 1 for key in posvar), 'NoVar')

            maxdual.update()
            maxdual.optimize()

            if maxdual.Status == 3:
                return model

            totalsum = maxdual.ObjVal
            print(totalsum)
            for var in maxdual.getVars():
                if var.x > 0.001:
                    print(var.varName, var.x)

            if totalsum > 0.01:
                nam = '_var' + str(len(model.getVars())) + '_'

                # Make a new variable corresponding to the matching and round
                model.addVar(vtype="CONTINUOUS", lb=0, ub=1, name=nam)
                model.update()
                print(nam)

                # put the edges in the matching into a list
                es = []
                for var in maxdual.getVars():
                    if var.varName[0] != 'v':
                        if var.x > 0.5:
                            print(var)
                            es.append(edges[int(var.varName[2:-1])])

                matchcost = 0

                # Update the constraints and objective value
                con = model.getConstrs()
                print(r, con[r])
                # Add the variable to its round constraint
                model.chgCoeff(con[r], model.getVars()[-1], 1)
                # For every match in the matching, add the variable to the match constraint
                for edge in es:
                    matchcost += costdict.get((edge[0], edge[1], r), 0)
                    index = edges.index(edge) + n - 1
                    print(edge, con[index])
                    model.chgCoeff(con[index], model.getVars()[-1], 1)
                # Add the variable to the objective function
                model.getVars()[-1].Obj = matchcost
                # Solve the primal again
                model.update()
                model.optimize()
                break
    return model


def solve_inst(cg, n, model, edges, costdict):
    print('Solving')
    integer = 0
    it = 0
    objs = [model.getObjective().getValue()]
    while not integer:
        it += 1
        vardict, posvar = SmallCG.construct_vardict(model, edges)
        odd_cycle, weight = find_odd_cycle(vardict)
        print('found cycle', odd_cycle)
        if weight > 0.999:
            return False, it, objs
        model.addConstr(
            (quicksum(model.getVarByName("_var" + var + "_") for var in odd_cycle) <= (len(odd_cycle) - 1) / 2), 'c')
        model.update()
        model.optimize()
        if (cg):
            model = solve_dual(model, n, edges, costdict, posvar)
        objs.append(model.ObjVal)
        for var in model.getVars():
            if var.x > 0.001:
                print(var)
        integer = 1
        for var in model.getVars():
            if 0.001 < var.x < 0.999:
                integer = 0
    return True, it, objs


def solve_indep_set(cg=False):
    # list with models that have an initially fractional solution
    number_viol = []
    # list with models that have an initially integral solution
    int_sol = []
    # list with models that have an initially fractional solution, but no violated inequalities
    no_found = []
    # dictionary that keeps track of the number of violated inequalities and iterations of a fractional model
    nrof_viol = dict()

    data_per_n = dict()

    # Iterate over all possible test cases
    # for s in ['06', '12', '18']:
    for s in ['06', '12']:
        G = nx.complete_graph(int(s))
        matchings = MakeInstance.make_matchings(G, [], [])
        for p in range(5, 10):
            nrof_nonint_inst = 0
            nrof_solved = 0
            avg_rel_obj = 0
            avg_cg_obj = 0
            algtime = 0
            for k in range(50):
                print("DIT IS K, S, P", k, s, p)
                if k < 10:
                    file = 'bin0' + s + '_0' + str(p) + '0_00' + str(k) + '.srr.lp'
                else:
                    file = 'bin0' + s + '_0' + str(p) + '0_0' + str(k) + '.srr.lp'
                if not cg:
                    model = MakeInstance.get_model(list(G.edges), matchings, file, s)
                    model.write('mod.lp')
                else:
                    model = read('results_rootrelaxation/' + file)
                model.Params.LogToConsole = 0
                for var in model.getVars():
                    if 'var' not in var.VarName:
                        model.remove(var)
                model.update()
                model.optimize()
                obj = model.getObjective().getValue()
                for var in model.getVars():
                    if var.x > 0:
                        print(var)
                # Get the parameters and dictionaries
                n, edges = SmallCG.get_params(model)
                vardict, posvar = SmallCG.construct_vardict(model, edges)
                costdict, n = SmallCG.make_costdict(file)

                for var in model.getVars():
                    if var.x > 0:
                        print(var.varName, vardict[var.varName])

                # If only n-1 variables are non-zero, there is an integer solution
                if len(vardict) == n - 1:
                    int_sol.append((s, p, k))
                else:
                    # Find independent set violations
                    odd_cycle, weight = find_odd_cycle(vardict)
                    # Try to solve the instance
                    starttime = time.time()
                    solved, it, objs = solve_inst(cg, n, model, edges, costdict)
                    algtime += time.time() - starttime
                    nrof_nonint_inst += 1
                    if solved:
                        nrof_solved += 1
                    avg_rel_obj += obj
                    avg_cg_obj += objs[-1]
                    if weight < 0.999:
                        number_viol.append((s, p, k))
                        nrof_viol[(s, p, k)] = {'nrof_iterations': it, 'solved': solved}
                    else:
                        no_found.append((s, p, k))
            algtime = algtime / nrof_nonint_inst
            avg_rel_obj = round(avg_rel_obj / nrof_nonint_inst, 3)
            avg_cg_obj = round(avg_cg_obj / nrof_nonint_inst, 3)
            data_per_n[(s, p)] = {'nrof fractional instances': nrof_nonint_inst,  'nrof solved': nrof_solved, 'avg obj rel' : avg_rel_obj, 'avg obj cg' : avg_cg_obj, 'avg time': algtime}
    print(len(number_viol), len(int_sol), no_found, nrof_viol, data_per_n)
    df = pd.DataFrame.from_dict(data_per_n, orient="index")
    print(df.to_latex())
    return len(number_viol), len(int_sol), no_found, nrof_viol, data_per_n

