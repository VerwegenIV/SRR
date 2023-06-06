from gurobipy import *
import math
import SmallCG
import networkx as nx
import MakeInstance
import pandas as pd


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


def solve_dual(model, n, edges, costdict, cycles):
    totalsum = 1

    while totalsum > 0.001:
        # dual vars for rounds constraints
        dual_r = model.Pi[:n - 1]
        # dual vars for match constraints
        dual_m = model.Pi[n - 1:n + int((n - 1) * n / 2) - 1]
        # dual vars for the cut constraints
        dual_c = model.Pi[n + int((n - 1) * n / 2) - 1:]

        for r in range(n - 1):
            print('Trying round', r)
            G = nx.Graph()
            for i in range(len(edges)):
                G.add_edge(edges[i][0], edges[i][1], weight=dual_m[i] - costdict.get((edges[i][0], edges[i][1], r), 0))
            match = nx.max_weight_matching(G, True)
            weight = 0
            for m in match:
                weight += G[m[0]][m[1]]['weight']
            print(match, weight)
        totalsum = 0


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
            model = solve_dual(model, n, edges, costdict, [])
        objs.append(model.getObjective().getValue())
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
    for s in ['06']:
        for p in range(5, 10):
            avg_iter = 0
            nrof_nonint_inst = 0
            avg_obj_incr = 0
            nrof_solved = 0
            for k in range(50):
                print("DIT IS K, S, P", k, s, p)
                if k < 10:
                    file = 'bin0' + s + '_0' + str(p) + '0_00' + str(k) + '.srr.lp'
                else:
                    file = 'bin0' + s + '_0' + str(p) + '0_0' + str(k) + '.srr.lp'
                if not cg:
                    model = MakeInstance.get_model(file, s)
                    model.write('mod.lp')
                else:
                    model = read('results_rootrelaxation/' + file)
                model.Params.LogToConsole = 0
                for var in model.getVars():
                    if 'var' not in var.VarName:
                        model.remove(var)
                model.update()
                model.optimize()
                for var in model.getVars():
                    if var.x > 0:
                        print(var)
                # Get the parameters and dictionaries
                n, edges = SmallCG.get_params(model)
                vardict, posvar = SmallCG.construct_vardict(model, edges)
                costdict = SmallCG.make_costdict(file)

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
                    solved, it, objs = solve_inst(cg, n, model, edges, costdict)
                    nrof_nonint_inst += 1
                    avg_iter += it
                    if solved:
                        nrof_solved += 1
                    avg_obj_incr += objs[-1] - objs[0]
                    if weight < 0.999:
                        number_viol.append((s, p, k))
                        nrof_viol[(s, p, k)] = {'nrof_iterations': it, 'solved': solved}
                    else:
                        no_found.append((s, p, k))
            avg_obj_incr = avg_obj_incr / nrof_nonint_inst
            avg_iter = avg_iter / nrof_nonint_inst
            data_per_n[(s, p)] = {'nrof fractional instances': nrof_nonint_inst, 'average objective increase':
                avg_obj_incr, 'average nrof iterations': avg_iter, 'nrof solved': nrof_solved}
    print(len(number_viol), len(int_sol), no_found, nrof_viol, data_per_n)
    df = pd.DataFrame.from_dict(data_per_n, orient="index")
    print(df.to_latex())
    return len(number_viol), len(int_sol), no_found, nrof_viol, data_per_n

