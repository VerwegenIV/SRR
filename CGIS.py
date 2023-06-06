from gurobipy import *
import math
import SmallCG
import IndepSet
import MakeInstance
import networkx as nx
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


def solve(n, model, edges):
    print('Solving')
    integer = 0
    it = 0
    nrof_cuts = 0
    objs = [model.getObjective().getValue()]
    while not integer:
        it += 1
        vardict, posvar = SmallCG.construct_vardict(model, edges)

        odd_cycle, weight = find_odd_cycle(vardict)
        if weight > 0.999:
            # Find violated inequalities
            alpha, beta, gamma = SmallCG.a_b_g(vardict)
            violated = SmallCG.find_violation(n, alpha, beta, gamma)
            if not violated:
                return False, it, objs

            # CG constraints
            for (mm1, mm2, r1) in violated:
                rsm = list(range(n - 1))
                rsm.remove(r1)
                model.addConstr(((quicksum(model.getVarByName(var) for var in posvar if (
                        mm1 not in posvar[var]['edges'] and mm2 not in posvar[var]['edges'] and r1 == posvar[var]['r']))
                                  - quicksum(model.getVarByName(var) for var in posvar if (
                                mm1 in posvar[var]['edges'] and mm2 in posvar[var]['edges']
                                and posvar[var]['r'] in rsm))) >= 0), 'cut' + str(nrof_cuts))
                nrof_cuts += 1
        else:
            # cycle constraints
            model.addConstr(
                (quicksum(model.getVarByName("_var" + var + "_") for var in odd_cycle) <= (len(odd_cycle) - 1) / 2),
                'c')


        # # Find violated inequalities
        # alpha, beta, gamma = SmallCG.a_b_g(vardict)
        # violated = SmallCG.find_violation(n, alpha, beta, gamma)
        # if not violated:
        #     odd_cycle, weight = find_odd_cycle(vardict)
        #     # cycle constraints
        #     model.addConstr(
        #         (quicksum(model.getVarByName("_var" + var + "_") for var in odd_cycle) <= (len(odd_cycle) - 1) / 2), 'c')
        #     if weight > 0.999:
        #         return False, it, objs
        # else:
        #     # CG constraints
        #     for (mm1, mm2, r1) in violated:
        #         rsm = list(range(n - 1))
        #         rsm.remove(r1)
        #         model.addConstr(((quicksum(model.getVarByName(var) for var in posvar if (
        #                 mm1 not in posvar[var]['edges'] and mm2 not in posvar[var]['edges'] and r1 == posvar[var]['r']))
        #                           - quicksum(model.getVarByName(var) for var in posvar if (
        #                         mm1 in posvar[var]['edges'] and mm2 in posvar[var]['edges']
        #                         and posvar[var]['r'] in rsm))) >= 0), 'cut' + str(nrof_cuts))
        #         nrof_cuts += 1
        model.update()
        model.optimize()
        objs.append(model.getObjective().getValue())
        integer = 1
        for var in model.getVars():
            if var.x > 0.001:
                print(var)
                if var.x < 0.999:
                    integer = 0
    return True, it, objs


def solve_cgis():
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
                model = MakeInstance.get_model(file, s)
                model.write('mod.lp')
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
                    solved, it, objs = solve(n, model, edges)
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

