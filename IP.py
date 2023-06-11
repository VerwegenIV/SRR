from gurobipy import *
import pandas as pd


def solve_IP():
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
    for s in ['06', '12']:
        for p in range(5, 10):
            nrof_nonint_inst = 0
            avg_int_obj = 0
            avg_rel_obj = 0
            for k in range(50):
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
                obj = model.ObjVal

                nonint = False
                for var in model.getVars():
                    if 0.0001 < var.x < 0.99999:
                        nonint = True
                    var.setAttr(GRB.Attr.VType, GRB.BINARY)

                model.update()
                model.optimize()
                if model.Status == 3:
                    continue

                print(model.Status)

                if nonint:
                    nrof_nonint_inst += 1
                    avg_int_obj += model.ObjVal
                    avg_rel_obj += obj
            avg_int_obj = avg_int_obj / round(nrof_nonint_inst, 2)
            avg_rel_obj = avg_rel_obj / round(nrof_nonint_inst, 2)
            data_per_n[(s, p)] = {'nrof fractional instances': nrof_nonint_inst, 'avg int obj': avg_int_obj}
    print(len(number_viol), len(int_sol), no_found, nrof_viol, data_per_n)
    df = pd.DataFrame.from_dict(data_per_n, orient="index")
    print(df.to_latex())
    return len(number_viol), len(int_sol), no_found, nrof_viol, data_per_n
