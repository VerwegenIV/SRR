from gurobipy import *
import networkx as nx


def make_matchings(G, match, matchings):
    if not G.edges:
        matchingset = set(frozenset(i) for i in matchings)
        if set(match) not in matchingset:
            matchings.append(match)
    else:
        for edge in G.edges:
            H = nx.Graph(G)
            m = match.copy()
            m.append(edge)
            H.remove_nodes_from(edge)
            make_matchings(H, m, matchings)
    return matchings


def parse_data(data):
    edge_one = []
    for d in data:
        l = re.split(' ', d)
        e = []
        for num in range(len(l)):
            if l[num]:
                e.append(l[num])
        if int(e[0]) < int(e[1]):
            edge_one.append((int(e[0]), int(e[1]), int(e[2])))
        else:
            edge_one.append((int(e[1]), int(e[0]), int(e[2])))
    edge_one = list(set(edge_one))
    return edge_one


def get_vars(costs, n):
    matchcost = dict()
    G = nx.complete_graph(n)
    matchings = make_matchings(G, [], [])
    for r in range(n - 1):
        for i in range(len(matchings)):
            for j in range(len(matchings[i])):
                if (sorted(matchings[i][j])[0], sorted(matchings[i][j])[1], r) in costs:
                    matchcost[r * len(matchings) + i] = matchcost.get(r * len(matchings) + i, 0) + 1
                else:
                    matchcost[r * len(matchings) + i] = matchcost.get(r * len(matchings) + i, 0)
    return list(G.edges), matchings, matchcost


def make_model(n, edges, matchings, costs):
    model = Model('SRR')
    model.ModelSense = GRB.MINIMIZE

    y = {}
    for key in costs:
        y[key] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="_var%d_" % key, obj=costs[key])

    for r in range(n - 1):
        model.addConstr((quicksum(y[r * len(matchings) + i] for i in range(len(matchings))) == 1), 'round[%d]' % r)

    for e in range(len(edges)):
        model.addConstr((quicksum(y[r*len(matchings) + i] for r in range(n-1) for i in range(len(matchings)) if edges[e] in
                                  matchings[i]) == 1), 'match[%d' % e)
    model.write('newmod.lp')
    return model


def get_model(file, s):
    n = int(s)
    data = open('Instances/' + file[:-3])
    data = data.read()
    data = data[len(str(n)) + 1:-1].split('\n')
    costs = parse_data(data)
    edges, matchings, costs = get_vars(costs, n)
    model = make_model(n, edges, matchings, costs)
    return model
