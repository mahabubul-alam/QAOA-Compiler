import json
import networkx as nx

n = 4
p = 0.8
seed = 0
g = nx.erdos_renyi_graph(n,p,seed=seed)

dic = {}

for i, edge in enumerate(g.edges()):
    dic['{}'.format(i+1)] = "{}".format(edge)

print(dic)

with open('QAOA_circ.json','w') as fp:
    json.dump(dic,fp)

