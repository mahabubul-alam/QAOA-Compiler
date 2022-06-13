"""
This script can be used/modified to create
custom quantum hardware configuration file
(an input to the compiler).
"""

import json
import random
import networkx as nx

#List of the single-qubit native gates supported by the hardware
single_q_gates = ['id','sx','x','u3']
#Two qubit native gate supported by the hardware,
#current compiler implementation assumes that the hardware
#supports only 1 two-qubit operation (cx)
two_q_gate = ['cx']

Q = 20
NEIGHBOR = 4
SEED = 0

coupling_graph = nx.random_regular_graph(NEIGHBOR,Q,seed=SEED)
dic = {'1Q': single_q_gates, '2Q': two_q_gate}
for gate in single_q_gates:
    tdic = {f'{q}': 1 for q in range(Q)}
    dic[f'{gate}'] = tdic

tdic = {
    f'({edge[0]},{edge[1]})': random.uniform(0.96, 0.99)
    for edge in coupling_graph.edges()
}

dic[f'{two_q_gate[0]}'] = tdic

with open('QC.json','w') as fp:
    json.dump(dic,fp,indent=4)
