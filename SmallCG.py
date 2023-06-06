from gurobipy import *
import math
import pandas as pd
import MakeInstance


def make_costdict(file):
    # Pick an instance
    data = open('Instances/' + file[:-3])
    # Since the instant is read as a string, we need to manually interpret it
    data = data.read()
    data = data.split('\n')
    # First entry is nrof nodes
    n = int(data[0])
    data = data[1:-1]
    costdict = dict()
    # Every entry in data is a string that contains the endpoints of an edge, a round number and the value 1.
    for d in data:
        # split them on spaces
        entries = re.split(' ', d)
        # make an edge entry with e[0] = left endpoint; e[1] = right endpoint; e[2] = round.
        e = []
        # for every non space entry
        for num in range(len(entries)):
            if entries[num]:
                e.append(entries[num])
        # Make sure that the left endpoint is always smaller
        if int(e[0]) < int(e[1]):
            costdict[(int(e[0]), int(e[1]), int(e[2]))] = 1
        else:
            costdict[(int(e[1]), int(e[0]), int(e[2]))] = 1
    print('COST', type(costdict))
    return costdict, n


# Get the parameters of the loaded model
def get_params(model):
    n = int(math.sqrt(2 * len(model.getConstrs()) + 2.25) - 0.5)
    edges = []
    for i in range(n - 1):
        for j in range(i + 1, n):
            edges.append((i, j))
    return n, edges


def construct_vardict(model, edges):
    # Initialise two dicts: one for every variable and one for positive variables
    vardict = dict()
    pos_vardict = dict()
    for var in model.getVars():
        # The matching that var represents
        es = []
        # The column of the variable; i.e. in which constraints its included
        c = model.getCol(var)
        # For a loaded in model, some non-variables are recognised as variables, this if-statement filters those out
        if c.size() <= 1:
            continue
        r = 0
        # For every constraint the variable is in
        for i in range(c.size()):
            name = c.getConstr(i).ConstrName
            # if its a round constraint, get the round
            if 'round' in name:
                r = int(name[6:-1])
            # if its a match constraint, get the edge
            elif 'match' in name:
                es.append(edges[int(name[6:])])
        # If the value is fractional, add it to the positive variables
        if var.x > 0:
            pos_vardict[var.VarName] = {'r': r, 'edges': es, 'value': var.x}
        # Add the variable with its value, the round and the matching it represents to the vardict.
        vardict[var.VarName] = {'r': r, 'edges': es, 'value': var.x}
    # If there are fractional variables, return it. Otherwise return only vardict.
    if pos_vardict:
        return pos_vardict, vardict
    else:
        return False, vardict


# Define the alpha-beta-gamma dictionaries
def a_b_g(vardict):
    alpha = dict()
    beta = dict()
    gamma = dict()
    # For every variable in the vardict
    for v in vardict:
        var = vardict[v]
        mat = vardict[v]['edges']
        # for every match in the matching
        for i in range(len(mat)):
            # Add the value of the variable y_M,r to alpha[m,r]
            if (mat[i], var['r']) not in alpha:
                alpha[(mat[i], var['r'])] = var['value']
            else:
                alpha[(mat[i], var['r'])] += var['value']
            # for every combination of matches in the matching
            for j in range(i + 1, len(mat)):
                # Add the value of the variable y_M,r to beta[m,m']
                if (mat[i], mat[j]) not in beta:
                    beta[(mat[i], mat[j])] = var['value']
                else:
                    beta[(mat[i], mat[j])] += var['value']
                # Add the value of the variable y_M,r to gamma[m,m',r]
                if (mat[i], mat[j], var['r']) not in gamma:
                    gamma[(mat[i], mat[j], var['r'])] = var['value']
                else:
                    gamma[(mat[i], mat[j], var['r'])] += var['value']
    return alpha, beta, gamma


# Find all violated inequalities using the alpha-beta-gamma dictionaries.
def find_violation(n, alpha, beta, gamma):
    violated = []
    for (m1, m2) in beta:
        for r in range(n - 1):
            if alpha.get((m1, r), 0) + alpha.get((m2, r), 0) + beta[(m1, m2)] - 2 * gamma.get((m1, m2, r), 0) > 1.001:
                violated.append((m1, m2, r))
    return violated


def solve_dual(model, n, edges, costdict, violated):
    # Initialise totalsum
    totalsum = 1
    print(costdict)
    print(violated)

    # While we find variables with negative reduced cost
    while totalsum > 0.01:
        # Create empty dictionaries to store the coefficients of (pairs of) variables
        dual_edge = dict()
        dual_pair = dict()

        # dual vars for rounds constraints
        dual_r = model.Pi[:n - 1]
        # dual vars for match constraints
        dual_m = model.Pi[n - 1:n + int((n - 1) * n / 2) - 1]
        # dual vars for the cut constraints
        dual_c = model.Pi[n + int((n - 1) * n / 2) - 1:]

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

            maxdual.update()
            maxdual.optimize()

            totalsum = maxdual.objVal
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
                # More debug code
                model.write('test.lp')
                totalsum = 0
                break
    return model


# Solve the models with column generation, adding all violated inequalities every iteration.
def solve_mod(cg, model, n, edges, costdict):
    vardict, posvar = construct_vardict(model, edges)
    # Find violated inequalities
    alpha, beta, gamma = a_b_g(vardict)
    violated = find_violation(n, alpha, beta, gamma)
    all_violated = violated.copy()
    # iteration counter
    it = 0
    nrof_cuts = 0
    # keep track of objective values through the run
    objs = [model.getObjective().getValue()]

    # While there are violated inequalities
    while violated:
        for var in model.getVars():
            if var.x > 0.001:
                print('VAR', var, vardict[var.varName])
        print('VIOLATED', violated)
        # Add every violated inequality to the model and solve it again
        for (mm1, mm2, r1) in violated:
            rsm = list(range(n - 1))
            rsm.remove(r1)
            model.addConstr(((quicksum(model.getVarByName(var) for var in posvar if (
                    mm1 not in posvar[var]['edges'] and mm2 not in posvar[var]['edges'] and r1 == posvar[var]['r']))
                              - quicksum(model.getVarByName(var) for var in posvar if (
                            mm1 in posvar[var]['edges'] and mm2 in posvar[var]['edges']
                            and posvar[var]['r'] in rsm))) >= 0), 'cut' + str(nrof_cuts))
            nrof_cuts += 1
        model.update()
        model.optimize()
        if (cg):
            model = solve_dual(model, n, edges, costdict, all_violated)
        objs.append(model.getObjective().getValue())

        # Update the variables in the vardict.
        vardict, posvar = construct_vardict(model, edges)
        if len(vardict) == n - 1:
            break
        # Find violated inequalities
        alpha, beta, gamma = a_b_g(vardict)
        violated = find_violation(n, alpha, beta, gamma)
        all_violated.extend(violated)
        it += 1

    # If there are no more violated inequalities return the number of iterations and whether the model is solved.
    for var in model.getVars():
        if 0.001 < var.x < 0.999:
            return it, objs, False
    return it, objs, True


def solve_instances(cg=False):
    # list with models that have an initially fractional solution
    number_viol = []
    # list with models that have an initially integer solution
    int_sol = []
    # list with models that have an initially fractional solution, but no violated inequalities
    no_found = []
    # dictionary that keeps track of the number of violated inequalities and iterations of a fractional model
    nrof_viol = dict()
    # dictionary that keeps the average iterations, average objective increase, and nrof instances solved per n
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
                print("DIT IS K,S,P", k, s, p)
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
                # Get the parameters and dictionaries
                n, edges = get_params(model)
                costdict, np = make_costdict(file)
                print(type(costdict))
                vardict, posvar = construct_vardict(model, edges)
                for var in model.getVars():
                    if var.x > 0.0001:
                        print(var, vardict[var.varName], var.Obj)
                # If only n-1 variables are non-zero, there is an integer solution
                if len(vardict) == n - 1:
                    int_sol.append((s, p, k))
                else:
                    # Solve the problems with integer solution to optimality
                    nrof_iter, objs, solved = solve_mod(cg, model, n, edges, costdict)
                    nrof_nonint_inst += 1
                    avg_iter += nrof_iter
                    if solved:
                        nrof_solved += 1
                    avg_obj_incr += objs[-1] - objs[0]
                    alpha, beta, gamma = a_b_g(vardict)
                    violated = find_violation(n, alpha, beta, gamma)
                    if violated:
                        number_viol.append((s, p, k))
                        nrof_viol[(s, p, k)] = {'nrof_viol': len(violated), 'nrof_iter': nrof_iter, 'objectives': objs, 'solved': solved}
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
