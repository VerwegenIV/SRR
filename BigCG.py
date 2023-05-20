from gurobipy import *
import math
import itertools
import SmallCG


def find_subsets(s, p, n, mts):
    allsets = set(itertools.combinations(s, p))
    # Remove subsets that are not disjoint.
    allsets2 = allsets.copy()
    if mts:
        for mat in allsets:
            for i in range(n):
                count = 0
                for e in mat:
                    if i in e:
                        count += 1
                if count > 1:
                    allsets2.remove(mat)
                    break
    return allsets2


def frac_edge_round(vardict):
    # Determine the edges and rounds that are fractional
    frac_edges = set()
    frac_round = set()
    for v in vardict:
        es = vardict[v]['edges']
        frac_round.add(vardict[v]['r'])
        for edge in es:
            frac_edges.add(edge)
    return frac_edges, frac_round


def solve_inst(model):
    n, edges = SmallCG.get_params(model)
    vardict, posvar = SmallCG.construct_vardict(model, edges)
    frac_edges, frac_round = frac_edge_round(vardict)
    viol = []
    for i in range(3, int(n / 2) + 1):
        # Find all possible combinations of i disjoint matches
        matches = find_subsets(frac_edges, i, n, True)
        if i - 1 > len(frac_round):
            if (len(frac_round) + i) % 2 == 0:
                ub = len(frac_round) - 1
            else:
                ub = len(frac_round)
        else:
            ub = i - 1
        for mats in matches:
            # find all possible combinations of j rounds
            for j in range(ub, 0, -2):
                rounds = find_subsets(frac_round, j, n, False)
                # define the right-hand-side
                rhs = math.ceil((i + j) / 2)
                # for every subset of fractional rounds
                for rnds in rounds:
                    totalsum = 0
                    for v in posvar:
                        var = posvar[v]
                        # Sum up the left-hand-side
                        if var['r'] in rnds:
                            coeff = math.ceil((1 + len(set(var['edges']).intersection(mats))) / 2)
                            totalsum += coeff * var['value']
                        else:
                            coeff = math.ceil((len(set(var['edges']).intersection(mats))) / 2)
                            totalsum += coeff * var['value']
                    # If the inequality is violated, add it to the list
                    if totalsum < rhs:
                        viol.append((mats, rnds))
    return viol

def solve_big_cg():
    nr_viol, int_sol, no_found, violations = SmallCG.solve_instances()
    viol_inst = dict()
    for inst in no_found:
        if inst[2] < 10:
            file = 'bin0' + inst[0] + '_0' + str(inst[1]) + '0_00' + str(inst[2]) + '.srr.lp'
        else:
            file = 'bin0' + inst[0] + '_0' + str(inst[1]) + '0_0' + str(inst[2]) + '.srr.lp'
        model = read('results_rootrelaxation/' + file)
        model.update()
        model.optimize()
        viol_inst[inst] = solve_inst(model)
    return viol_inst
