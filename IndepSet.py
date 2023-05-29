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
        vars_in_0.append(dict_keys[-1][4:-1])
    if len(vars_in_0) < 2:
        return odd_cycles
    for node in vars_in_0:
        v_in_0 = vars_in_0.copy()
        v_in_0.remove(node)
        H = nx.Graph(G)
        H.remove_nodes_from(v_in_0)
        print(H.edges)
        cycles = cycle_sub(n, node, 'a', [node], H, set())
        odd_cycles.extend(cycles)
    return odd_cycles


def cycle_sub(n, node, prev_node, current_c, G, cycles):
    print(node, prev_node, current_c)
    if len(current_c) == n - 1 and current_c[0] in list(G.adj[node]):
        cycles.add(frozenset(current_c))
        return cycles
    c_c = current_c.copy()
    if len(c_c) > 1:
        c_c.remove(prev_node)
    if len(set(c_c).intersection(set(G.adj[node]))) > 0:
        return cycles
    adj_list = list(G.adj[node])
    if prev_node != 'a':
        adj_list.remove(prev_node)
    for new_node in adj_list:
        H = nx.Graph(G)
        current_c.append(new_node)
        for other_node in list(G.adj[new_node]):
            if G.nodes[other_node]['r'] == G.nodes[new_node]['r']:
                H.remove_node(other_node)
        cycles = cycle_sub(n, new_node, node, current_c, H, cycles)
        current_c.remove(new_node)
    return cycles


def solve_inst(model, n, edges):
    print('Solving')
    old_triangles = []
    old_cycles = set()
    integer = 0
    while not integer:
        vardict, posvar = SmallCG.construct_vardict(model, edges)
        triangles = find_triangles(vardict)
        odd_cycles = find_odd_cycle(n, vardict)
        if triangles == old_triangles or old_cycles == odd_cycles:
            return False
        old_triangles = triangles.copy()
        old_cycles = odd_cycles.copy()
        print('found triangles', triangles)
        print('found cycles', odd_cycles)
        if not triangles and not odd_cycles:
            return False
        for triangle in triangles:
            for var in triangle:
                print(var, model.getVarByName(var))
            model.addConstr((quicksum(model.getVarByName(var) for var in triangle) <= 1), 't')
        for cycle in odd_cycles:
            model.addConstr((quicksum(model.getVarByName("_var" + var + "_") for var in cycle) <= (n - 2) / 2), 'c')
        model.update()
        model.optimize()
        for var in model.getVars():
            if var.x > 0.001:
                print(var)
        integer = 1
        for var in model.getVars():
            if 0.001 < var.x < 0.999:
                integer = 0
    return True


def solve_indep_set():
    # list with models that have an initially fractional solution
    number_viol = []
    # list with models that have an initially integral solution
    int_sol = []
    # list with models that have an initially fractional solution, but no violated inequalities
    no_found = []
    # dictionary that keeps track of the number of violated inequalities and iterations of a fractional model
    nrof_viol = dict()

    # Iterate over all possible test cases
    # for s in ['06', '12', '18']:
    for s in ['06']:
        for p in range(5, 10):
            for k in range(50):
                print("DIT IS K, S, P", k, s, p)
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
                    triangles = find_triangles(vardict)
                    odd_cycles = find_odd_cycle(n, vardict)
                    # Try to solve the instance
                    solved = solve_inst(model, n, edges)
                    if triangles or odd_cycles:
                        number_viol.append((s, p, k))
                        nrof_viol[(s, p, k)] = {'nrof_triangles': len(triangles), 'nrof_cycles': len(odd_cycles),
                                                'solved': solved}
                    else:
                        no_found.append((s, p, k))
    return len(number_viol), len(int_sol), no_found, nrof_viol
