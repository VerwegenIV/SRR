from gurobipy import *
import math
import SmallCG
import networkx as nx


def find_triangles(vardict):
    dict_keys = list(vardict.keys())
    triangles = []
    for i in range(len(dict_keys) - 2):
        var1 = vardict[dict_keys[i]]
        for j in range(i + 1, len(dict_keys) - 1):
            var2 = vardict[dict_keys[j]]
            for k in range(j + 1, len(dict_keys)):
                var3 = vardict[dict_keys[k]]
                if len(set(var1['edges']).intersection(set(var2['edges']))) > 0 or var1['r'] == var2['r']:
                    if len(set(var1['edges']).intersection(set(var3['edges']))) > 0 or var1['r'] == var3['r']:
                        if len(set(var2['edges']).intersection(set(var3['edges']))) > 0 or var2['r'] == var3['r']:
                            if var1['value'] + var2['value'] + var3['value'] > 1:
                                triangles.append((dict_keys[i], dict_keys[j], dict_keys[k]))
    return triangles


def find_odd_cycle(n, vardict):
    dict_keys = list(vardict.keys())
    vars_in_0 = []
    odd_cycles = []
    G = nx.Graph()
    for var in vardict:
        G.add_node(var[4:-1], r=vardict[var]['r'])
    for i in range(len(dict_keys) - 1):
        var1 = dict_keys[i]
        if vardict[var1]['r'] == 0:
            vars_in_0.append(var1[4:-1])
        for j in range(i + 1, len(dict_keys)):
            var2 = dict_keys[j]
            if len(set(vardict[var1]['edges']).intersection(set(vardict[var2]['edges']))) > 0 or vardict[var1]['r'] == \
                    vardict[var2]['r']:
                G.add_edge(var1[4:-1], var2[4:-1])
    if vardict[dict_keys[-1]]['r'] == 0:
        vars_in_0.append(dict_keys[-1])
    for node in vars_in_0:
        v_in_0 = vars_in_0.copy()
        v_in_0.remove(node)
        H = nx.Graph(G)
        H.remove_nodes_from(v_in_0)
        cycles = cycle_sub(n, node, 'a', [node], H, {})
        odd_cycles.extend(cycles)
    print(odd_cycles)


def cycle_sub(n, node, prev_node, current_c, G, cycles):
    if len(current_c) == n - 1 and current_c[0] in list(G.adj[prev_node]):
        cycles.add(current_c)
        return cycles
    c_c = current_c.copy()
    if len(c_c) > 1:
        c_c.remove(prev_node)
    if len(set(c_c).intersection(set(G.adj[node]))) > 0:
        return cycles
    for new_node in list(G.adj[node]):
        H = nx.Graph(G)
        current_c.append(new_node)
        for other_node in list(G.adj[new_node]):
            if G.nodes[other_node]['r'] == G.nodes[new_node]['r']:
                H.remove_node(other_node)
        cycles = cycle_sub(n, new_node, node, current_c, H, cycles)
        current_c.remove(new_node)
    return cycles


def solve_indep_set():
    # list with models that have an initially fractional solution
    number_viol = []
    # list with models that have an initially integer solution
    int_sol = []
    # list with models that have an initially fractional solution, but no violated inequalities
    no_found = []
    # dictionary that keeps track of the number of violated inequalities and iterations of a fractional model
    nrof_viol = dict()

    # Iterate over all possible test cases
    # for s in ['06', '12', '18']:
    for s in ['06']:
        for p in range(5, 6):
            for k in range(32,33):
                print("DIT IS K,S,P", k, s, p)
                if k < 10:
                    file = 'bin0' + s + '_0' + str(p) + '0_00' + str(k) + '.srr.lp'
                else:
                    file = 'bin0' + s + '_0' + str(p) + '0_0' + str(k) + '.srr.lp'
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

                for var in model.getVars():
                    if var.x > 0:
                        print(var.varName, vardict[var.varName])

                # If only n-1 variables are non-zero, there is an integer solution
                if len(vardict) == n - 1:
                    int_sol.append((s, p, k))
                else:
                    # Find small independent set violations
                    # find_triangles(vardict)
                    find_odd_cycle(n, vardict)
                    # if violated:
                    #     number_viol.append((s, p, k))
                    #     nrof_viol[(s, p, k)] = {'nrof_viol': len(violated), 'nrof_iter': nrof_iter, 'objectives': objs}
                    # else:
                    #     no_found.append((s, p, k))
    # return len(number_viol), len(int_sol), no_found, nrof_viol
