from gurobipy import *
import math


def make_costdict(file):
    # Pick an instance
    data = open('Instances/'+file[:-3])
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
        l = re.split(' ', d)
        # make an edge entry with e[0] = left endpoint; e[1] = right endpoint; e[2] = round.
        e = []
        # for every non space entry
        for num in range(len(l)):
            if l[num]:
                e.append(l[num])
        # Make sure that the left endpoint is always smaller
        if int(e[0]) < int(e[1]):
            costdict[(int(e[0]), int(e[1]), int(e[2]))] = 1
        else:
            costdict[(int(e[1]), int(e[0]), int(e[2]))] = 1
    for i in range(n - 1):
        for j in range(i + 1, n):
            for r in range(n - 1):
                costdict[(i, j, r)] = costdict.get((i, j, r), 0)

    return costdict

# Get the parameters of the loaded model
def get_params(model):
    n = int(math.sqrt(2*len(model.getConstrs())+2.25) - 0.5)
    edges = []
    for i in range(n-1):
        for j in range(i+1,n):
            edges.append((i,j))
    return n,edges

def construct_vardict(model,edges):
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
            pos_vardict[var.VarName] = {'r':r, 'edges':es, 'value':var.x}
        # Add the variable with its value, the round and the matching it represents to the vardict.
        vardict[var.VarName] = {'r':r, 'edges':es, 'value':var.x}
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
            if (mat[i],var['r']) not in alpha:
                alpha[(mat[i],var['r'])] = var['value']
            else:
                alpha[(mat[i],var['r'])] += var['value']
            # for every combination of matches in the matching
            for j in range(i+1,len(mat)):
                # Add the value of the variable y_M,r to beta[m,m']
                if (mat[i],mat[j]) not in beta:
                    beta[(mat[i],mat[j])] = var['value']
                else:
                    beta[(mat[i],mat[j])] += var['value']
                # Add the value of the variable y_M,r to gamma[m,m',r]
                if (mat[i],mat[j],var['r']) not in gamma:
                    gamma[(mat[i],mat[j],var['r'])] = var['value']
                else:
                    gamma[(mat[i],mat[j],var['r'])] += var['value']
    return alpha, beta, gamma

# Find all violated inequalities using the alpha-beta-gamma dictionaries.
def find_violation(n, alpha,beta,gamma):
    violated = []
    for (m1,m2) in beta:
        for r in range(n-1):
            if alpha.get((m1,r),0) + alpha.get((m2,r),0) + beta[(m1,m2)] - 2*gamma.get((m1,m2,r),0) > 1.001:
                violated.append((m1,m2,r))
    return violated


def solve_dual2(model, n, edges, costdict, violated):
    print('SOLVING DUAL')
    # Store the maximum for every round in a dicht
    maxr = dict()
    # Initialise totalsum
    totalsum = 1

    while totalsum > 0.01:
        # dual vars for rounds constraints
        dual_r = model.Pi[:n - 1]
        # dual vars for match constraints
        dual_m = model.Pi[n - 1:n + int((n - 1) * n / 2) - 1]
        # dual vars for the cut constraints
        dual_c = model.Pi[n + int((n - 1) * n / 2) - 1:]
        for r in range(n - 1):
            print('Trying round', r)
            maxdual = Model("dual")
            maxdual.ModelSense = GRB.MAXIMIZE
            maxdual.Params.LogToConsole = 0
            x = maxdual.addVars(edges, vtype=GRB.BINARY, name='x')
            x_p = maxdual.addVars(range(len(violated)), vtype=GRB.INTEGER, name='x_p')
            x_n = maxdual.addVars(range(len(violated)), vtype=GRB.INTEGER, name='x_n')

            maxdual.setObjective(dual_r[r] + quicksum(
                (dual_m[edges.index(e)] - costdict.get((e[0], e[1], r), 0)) * x[e] for e in edges) + quicksum(
                dual_c[i] * x_p[i] for i in range(len(violated)) if r == violated[i][2]) - quicksum(
                dual_c[i] * x_n[i] for i in range(len(violated)) if r != violated[i][2]))
            maxdual.addConstrs(((quicksum(x[e] for e in edges if v in e) == 1) for v in range(n)), 'matching')
            maxdual.addConstrs(
                ((x_p[i] <= (2 - x[violated[i][0]] - x[violated[i][1]]) / 2) for i in range(len(violated))), 'pUB')
            maxdual.addConstrs(
                ((x_p[i] + 1 >= (2 - x[violated[i][0]] - x[violated[i][1]]) / 2 + 0.001) for i in range(len(violated))),
                'pLB')
            maxdual.addConstrs(((x_n[i] <= (x[violated[i][0]] + x[violated[i][1]]) / 2) for i in range(len(violated))),
                               'nUB')
            maxdual.addConstrs(
                ((x_n[i] + 1 >= (x[violated[i][0]] + x[violated[i][1]]) / 2 + 0.001) for i in range(len(violated))),
                'nLB')

            maxdual.update()
            maxdual.optimize()

            totalsum = maxdual.objVal
            print(totalsum)

            if totalsum > 0.01:
                #                 if len(model.getVars()) > 75:
                #                     return False
                nam = '_var' + str(len(model.getVars())) + '_'

                # Make a new variable corresponding to the matching and round
                m = model.addVar(vtype="CONTINUOUS", lb=0, ub=1, name=nam)
                model.update()
                print(nam)

                es = []
                for var in maxdual.getVars():
                    if var.x > 0.5:
                        print(var)
                        if 'x[' in var.VarName:
                            es.append((int(var.varName[2:var.VarName.index(',')]),
                                       int(var.Varname[var.VarName.index(',') + 1:-1])))

                matchcost = 0

                # Update the constraints and objective value
                con = model.getConstrs()
                print(r, con[r])
                model.chgCoeff(con[r], model.getVars()[-1], 1)
                for edge in es:
                    matchcost += costdict.get((edge[0], edge[1], r), 0)
                    index = edges.index(edge) + n - 1
                    print(edge, con[index])
                    model.chgCoeff(con[index], model.getVars()[-1], 1)
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
                model.getVars()[-1].Obj = matchcost
                model.update()
                model.optimize()
                break
    return model


# Solve the models with column generation, adding all violated inequalities every iteration.
def solve_mod_cg(model, n, edges, vardict, costdict):
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
    while (violated):
        # Add every violated inequality to the model and solve it again
        for (mm1, mm2, r1) in violated:
            rsm = list(range(n - 1))
            rsm.remove(r1)
            model.addConstr(((quicksum(model.getVarByName(var) for var in vardict if (
                        mm1 not in vardict[var]['edges'] and mm2 not in vardict[var]['edges'] and r1 == vardict[var][
                    'r'])) - quicksum(model.getVarByName(var) for var in vardict if (
                        mm1 in vardict[var]['edges'] and mm2 in vardict[var]['edges'] and vardict[var][
                    'r'] in rsm))) >= 0), 'cut' + str(nrof_cuts))
            nrof_cuts += 1
        model.update()
        model.optimize()
        model = solve_dual2(model, n, edges, costdict, all_violated)
        objs.append(model.getObjective().getValue())

        # Update the variables in the vardict.
        vardict, posvar = construct_vardict(model, edges)
        if len(vardict) == n - 1:
            break
        # Find violated inequalities
        alpha, beta, gamma = a_b_g(vardict)
        violated = find_violation(n, alpha, beta, gamma)
        print(violated)
        print(vardict)
        all_violated.extend(violated)
        it += 1

    # If there are no more violated inequalities, return the number of iterations, or 0 if the solution is still fractional.
    for var in model.getVars():
        if var.x > 0.001 and var.x < 0.999:
            return 0, objs
    return it, objs

def solve_instances():
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
    for s in ['06', '12']:
        for p in range(5, 10):
            for k in range(50):
                print("DIT IS K,S,P", k, s, p)
                if k < 10:
                    file = 'bin0' + s + '_0' + str(p) + '0_00' + str(k) + '.srr.lp'
                else:
                    file = 'bin0' + s + '_0' + str(p) + '0_0' + str(k) + '.srr.lp'
                model = Model('RRT')
                model = read('results_rootrelaxation/'+file)
                model.Params.LogToConsole = 0
                for var in model.getVars():
                    if 'var' not in var.VarName:
                        model.remove(var)
                model.update()
                model.optimize()
                # Get the parameters and dictionaries
                n, edges = get_params(model)
                costdict = make_costdict(file)
                vardict, posvar = construct_vardict(model, edges)
                # If only n-1 variables are non-zero, there is an integer solution
                if len(vardict) == n - 1:
                    int_sol.append((s, p, k))
                else:
                    # Solve the problems with integer solution to optimality
                    nrof_iter, objs = solve_mod_cg(model, n, edges, posvar, costdict)
                    # nrof_iter, objs = solve_mod(model, n, posvar)
                    # nrof_iter = 1
                    # objs = model.getObjective().getValue()
                    alpha, beta, gamma = a_b_g(vardict)
                    violated = find_violation(n, alpha, beta, gamma)
                    if violated:
                        number_viol.append((s, p, k))
                        nrof_viol[(s, p, k)] = {'nrof_viol': len(violated), 'nrof_iter': nrof_iter, 'objectives': objs}
                    else:
                        no_found.append((s, p, k))
    print(len(number_viol))
    print(len(int_sol))
    print(no_found)
    print(nrof_viol)


