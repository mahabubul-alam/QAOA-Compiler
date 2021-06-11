import networkx as nx
import json
import random

single_q_gates = ['id','sx','x','u3'] #List of the single-qubit native gates supported by the hardware
two_q_gate = ['cx'] #Two qubit native gate supported by the hardware, current implementation assumes that the hardware supports only 1 two-qubit operation

qubit = 20
coupling_each_qubit = 4
seed = 0
coupling_graph = nx.random_regular_graph(coupling_each_qubit,qubit,seed=seed)

dic = {}

dic['1Q'] = single_q_gates
dic['2Q'] = two_q_gate

for gate in single_q_gates:
    tdic = {}
    for q in range(qubit):
        tdic['{}'.format(q)] = 1
    dic['{}'.format(gate)] = tdic

gate = two_q_gate[0]
tdic = {}
for edge in coupling_graph.edges():
    tdic['({},{})'.format(edge[0],edge[1])] = random.uniform(0.96,0.99)
dic['{}'.format(gate)] = tdic

with open('QC.json','w') as fp:
    json.dump(dic,fp)

