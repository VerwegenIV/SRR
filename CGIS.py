from gurobipy import *
import SmallCG
import IndepSet
import MakeInstance
import networkx as nx
import pandas as pd
import time


def solve_dual(model, n, edges, costdict, violated, posvar):
    # Initialise totalsum
    totalsum = 1

    # While we find variables with negative reduced cost
    while totalsum > 0.01:
        # Create empty dictionaries to store the coefficients of (pairs of) variables
        dual_edge = dict()
        dual_pair = dict()

        # dual vars for rounds constraints
        dual_r = model.Pi[:n - 1]
        # dual vars for match constraints
        dual_m = model.Pi[n - 1:n + int((n - 1) * n / 2) - 1]

        # dual vars for the CG cut constraints
        dual_c = []
        for i in range(n + int((n - 1) * n / 2) - 1, len(model.Pi)):
            if 'cut' in model.getConstrs()[i].ConstrName:
                dual_c.append(model.Pi[i])

        # Fill the coefficient dictionaries
        for r in range(n - 1):
            for k in range(len(violated)):
                if violated[k][2] != r:
                    dual_pair[(edges.index(violated[k][0]), edges.index(violated[k][1]), r)] = \
                        dual_pair.get((edges.index(violated[k][0]), edges.index(violated[k][1]), r), 0) - dual_c[k]
                else:
                    dual_edge[(edges.index(violated[k][0]), r)] = dual_edge.get((edges.index(violated[k][0]), r), 0) - \
                                                                  dual_c[k]
                    dual_pair[(edges.index(violated[k][0]), edges.index(violated[k][1]), r)] = \
                        dual_pair.get((edges.index(violated[k][0]), edges.index(violated[k][1]), r), 0) + dual_c[k]
            for e in range(len(edges)):
                dual_edge[(e, r)] = dual_edge.get((e, r), 0) + dual_m[e] - costdict.get((edges[e][0], edges[e][1], r),
                                                                                        0)
        print(dual_edge, dual_pair)

        # For every round, construct the dual and check whether the matching with maximum value violates a constraint
        for r in range(n - 1):
            print('Trying round', r)
            maxdual = Model("dual")
            maxdual.ModelSense = GRB.MAXIMIZE
            maxdual.Params.LogToConsole = 0

            x = maxdual.addVars(range(len(edges)), vtype=GRB.BINARY, name='x')
            v = maxdual.addVars(posvar, vtype=GRB.INTEGER, name='v')
            # Objective of dual_r[r] + dual_edge[e,r] * x_e + dual_pair[e,f,r] * x_e * x_f
            maxdual.setObjective(dual_r[r] - quicksum(dual_c[k] for k in range(len(violated)) if violated[k][2] != r)
                                 + quicksum(dual_edge[(e, r)] * x[e] for e in range(len(edges)))
                                 + quicksum(
                dual_pair[(edges.index(violated[k][0]), edges.index(violated[k][1]), r)] * x[
                    edges.index(violated[k][0])]
                * x[edges.index(violated[k][1])] for k in range(len(violated))))
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

            # If the matching found violates the dual constraint
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
                # Add the variable to its round constraint
                model.chgCoeff(con[r], model.getVars()[-1], 1)
                # For every match in the matching, add the variable to the match constraint
                for edge in es:
                    matchcost += costdict.get((edge[0], edge[1], r), 0)
                    index = edges.index(edge) + n - 1
                    print(edge, con[index])
                    model.chgCoeff(con[index], model.getVars()[-1], 1)
                # Add the variable to the violated CG cuts where applicable
                for k in range(len(violated)):
                    index = n + int((n - 1) * n / 2) + k - 1
                    (m3, m4, r1) = violated[k]
                    print(m3, m4, r1, con[index])
                    if len(set(es).intersection({m3, m4})) == 2 and r1 != r:
                        print("negative", model.getVars()[-1])
                        model.chgCoeff(con[index], model.getVars()[-1], -1)
                    elif len(set(es).intersection({m3, m4})) == 0 and r1 == r:
                        print("positive", model.getVarByName(nam))
                        model.chgCoeff(con[index], model.getVarByName(nam), 1)
                # Add the variable to the objective function
                model.getVars()[-1].Obj = matchcost
                # Solve the primal again
                model.update()
                model.optimize()
                break
    return model


def solve(cg, n, model, edges, costdict):
    print('Solving')
    integer = 0
    it = 0
    nrof_cuts = 0
    objs = [model.getObjective().getValue()]
    all_violated = []
    while not integer:
        it += 1
        vardict, posvar = SmallCG.construct_vardict(model, edges)

        odd_cycle, weight = IndepSet.find_odd_cycle(vardict)
        if weight > 0.999:
            # Find violated inequalities
            alpha, beta, gamma = SmallCG.a_b_g(vardict)
            violated = SmallCG.find_violation(n, alpha, beta, gamma)
            all_violated.extend(violated)
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
                'oc')


        # # Find violated inequalities
        # alpha, beta, gamma = SmallCG.a_b_g(vardict)
        # violated = SmallCG.find_violation(n, alpha, beta, gamma)
        # all_violated.extend(violated)
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
        if (cg):
            model = solve_dual(model, n, edges, costdict, all_violated, posvar)
        objs.append(model.getObjective().getValue())
        integer = 1
        for var in model.getVars():
            if var.x > 0.001:
                print(var)
                if var.x < 0.999:
                    integer = 0
    return True, it, objs


def solve_cgis(cg = False):
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
                model.write('mod.lp')
                model.Params.LogToConsole = 0
                for var in model.getVars():
                    if 'var' not in var.VarName:
                        model.remove(var)
                model.update()
                model.optimize()
                obj = model.ObjVal
                print(model.getConstrs())
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
                    odd_cycle, weight = IndepSet.find_odd_cycle(vardict)
                    # Try to solve the instance
                    starttime = time.time()
                    solved, it, objs = solve(cg, n, model, edges, costdict)
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

