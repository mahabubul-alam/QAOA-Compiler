################################################################################
##############      This library has been created by      ######################
##############            Md Mahabubul Alam               ######################
############## https://mahabubul-alam.github.io/Personal/ ######################
##############         Graduate Student (Ph.D.)           ######################
##############   Department of Electrical Engineering     ######################
##############      Pennsylvania State University         ######################
##############         University Park, PA, USA           ######################
##############               mxa890@psu.edu               ######################
################################################################################
import math
import re
import os
from qiskit import QuantumCircuit, transpile
import commentjson as json
import networkx as nx
from qiskit.circuit import Parameter

class CompileQAOAQiskit:
    def __init__(self, Circuit_json=None, \
        QC_json=None, Config_json=None, Out_Circuit_File_Name='QAOA.qasm'):
        """
        This method initializes necessary config variables.
        """
        self.load_config(Config_json)
        self.output_file_name = Out_Circuit_File_Name
        self.extract_qc_data(QC_json)
        self.layer_zz_assignments = {}
        self.zz_graph = self.qaoa_zz_graph(Circuit_json)
        self.initial_layout = list(range(len(self.zz_graph.nodes())))
        self.circuit = None
        self.sorted_ops = None
        self.cost = 10e10
        self.final_map = []

    def load_config(self,config_json=None):
        """
        This method loads the variables in the config json file.
        """
        with open(config_json) as f:
            self.Config = json.load(f)

        if 'Target_p' in self.Config.keys():
            self.Target_p = int(self.Config['Target_p'])
        else:
            self.Target_p = 1
        if 'Packing_Limit' in self.Config.keys():
            self.Packing_Limit = float(self.Config['Packing_Limit'])
        else:
            self.Packing_Limit = 10e10
        if 'Route_Method' in self.Config.keys():
            self.Route_Method = self.Config['Route_Method']
        else:
            self.Route_Method = 'sabre'
        if 'Trans_Seed' in self.Config.keys():
            self.Trans_Seed = int(self.Config['Trans_Seed'])
        else:
            self.Trans_Seed = 0
        if 'Opt_Level' in self.Config.keys():
            self.Opt_Level = int(self.Config['Opt_Level'])
        else:
            self.Opt_Level = 1

    def extract_qc_data(self,QC_file=None):
        """
        This method extracts hardware information from the QC json file.
        """
        with open(QC_file) as f:
            self.QC_DATA = json.load(f)

        try:
            self.Native_2Q = eval(self.QC_DATA['2Q'].strip('"'))
        except:
            self.Native_2Q = self.QC_DATA['2Q']
        try:
            self.Native_1Q = eval(self.QC_DATA['1Q'].strip('"'))
        except:
            self.Native_1Q = self.QC_DATA['1Q']

        self.Basis_Gates = self.Native_2Q + self.Native_1Q
        self.Coupling_Map = []
        for key in self.QC_DATA[str(self.Native_2Q[0])].keys():
            n1, n2 = eval(key)[0], eval(key)[1]
            if [n1,n2] not in self.Coupling_Map:
                self.Coupling_Map.append([n1,n2])
            if [n2,n1] not in self.Coupling_Map:
                self.Coupling_Map.append([n2,n1])
        self.calc_qq_distances()

    def calc_qq_distances(self):
        """
        This method calculates pairwise qubit-qubit distances using the floyd_warshall algorithm.
        """
        self.Unweighted_Undirected_Coupling_Graph = nx.Graph()
        self.Weighted_Undirected_Coupling_Graph = nx.Graph()

        for key, value in self.QC_DATA[str(self.Native_2Q[0])].items():
            n1, n2, sp = eval(key)[0], eval(key)[1], float(value)
            self.Unweighted_Undirected_Coupling_Graph.add_edge(n1,n2)
            self.Weighted_Undirected_Coupling_Graph.add_edge(n1,n2,weight=1/sp)

        self.QQ_Distances = nx.floyd_warshall(self.Unweighted_Undirected_Coupling_Graph)
        self.Noise_Aware_QQ_Distances = nx.floyd_warshall(self.Weighted_Undirected_Coupling_Graph)

    def set_iterc_target(self,Target='GC_2Q'):
        """
        This method can be used to set the target of iterative compilation.
        """
        self.IterC_Target = Target

    def set_incrc_var_awarenes(self,Variation_Aware=False):
        """
        This method can be used to set variation awareness in incremental compilation.
        """
        self.IncrC_VarAwareness = Variation_Aware

    def qaoa_zz_graph(self,circ_json=None):
        """
        This method is used to create the MaxCut graph from the json file.
        """
        with open(circ_json) as f:
            data = json.loads(f.read())
        zz_graph = nx.Graph()
        for key, val in data.items():
            nodes = eval(val)
            zz_graph.add_edge(nodes[0],nodes[1])
        return zz_graph

    def final_mapping_ic(self,qiskit_ckt_object):
        """
        This method finds the output logical-to-physical qubit mapping in a compiled circuit block.
        """
        qiskit_ckt_object.qasm(filename='tempfile')
        os.system('grep measure tempfile | awk \'{print $2, $4}\' > temp')
        qreg_creg_map = open('temp','r').readlines()

        FMAP = {}
        for line in qreg_creg_map:
            elements = line.split(' ')
            physical_qs = elements[0]
            logical_qs = str(elements[1]).split(';')[0]

            physical_q = re.search('\[.*\]',physical_qs).group()
            physical_q = re.sub('\[|\]','',physical_q)

            logical_q = re.search('\[.*\]',logical_qs).group()
            logical_q = re.sub('\[|\]','',logical_q)

            FMAP[logical_q[0:]] = int(physical_q[0:])

        final_map = []
        for i in range(qiskit_ckt_object.width()-qiskit_ckt_object.num_qubits):
            final_map.append(FMAP[str(i)])
        os.system('rm temp')
        os.system('rm tempfile')
        self.final_map = final_map

    def set_initial_layout(self,target_layout=None):
        """
        This method is used to set initial layout before starting compilation of any circuit block.
        """
        if target_layout:
            self.initial_layout = target_layout

    def sort_zz_by_QQ_distances(self,Unsorted_ops=None):
        """
        This method sorts the ZZ operations based on the control-target distances for the current mapping.
        """
        Sorted_ops = []
        Swap_distances_ops = {}

        for op in Unsorted_ops:
            _physical_q1 = self.initial_layout[op[0]]
            _physical_q2 = self.initial_layout[op[1]]
            if not self.IncrC_VarAwareness:
                swap_dist = self.QQ_Distances[_physical_q1][_physical_q2]
            else:
                swap_dist = self.Noise_Aware_QQ_Distances[_physical_q1][_physical_q2]
            Swap_distances_ops[op] = swap_dist

        for op in Unsorted_ops:
            if not Sorted_ops:
                Sorted_ops.append(op)
                continue
            i = 0
            for sop in Sorted_ops:
                if Swap_distances_ops[op] < Swap_distances_ops[sop]:
                    Sorted_ops.insert(i,op)
                    break
                i = i + 1
            if i == len(Sorted_ops):
                Sorted_ops.append(op)
        self.sorted_ops = Sorted_ops

    def construct_single_layer_ckt_IC(self,p):
        """
        This method constructs a single layer of the circuit in incremental compilation.
        """
        n = len(self.zz_graph.nodes())
        qc = QuantumCircuit(n, n)
        for zz in self.layer_zz_assignments['L0']:
            n1 = zz[0]
            n2 = zz[1]
            gamma = Parameter('g{}_{}_{}'.format(p,n1,n2))
            qc.cx(n1, n2)
            qc.rz(gamma, n2)
            qc.cx(n1, n2)
        qc.measure(range(n), range(n))
        trans_ckt = transpile(qc, coupling_map = self.Coupling_Map, \
                basis_gates = self.Basis_Gates, initial_layout = self.initial_layout, \
                optimization_level = self.Opt_Level, seed_transpiler = self.Trans_Seed, routing_method = self.Route_Method)
        self.circuit =  trans_ckt

    def incremental_compilation(self):
        """
        This method is used for incremental compilation.
        """
        logical_n = len(self.zz_graph.nodes())
        physical_n = len(self.Unweighted_Undirected_Coupling_Graph.nodes())
        IncrC_qc = QuantumCircuit(physical_n, logical_n)

        for i in range(logical_n):
            IncrC_qc.h(self.initial_layout[i])

        for p in range(1,self.Target_p+1):
            Remaining_ops = self.zz_graph.edges()
            while Remaining_ops:
                self.sort_zz_by_QQ_distances(Unsorted_ops=Remaining_ops)
                Sorted_ops = self.sorted_ops
                self.instruction_parallelization(Sorted_ops,Single_Layer=True)

                Remaining_ops = self.layer_zz_assignments['R']

                self.construct_single_layer_ckt_IC(p)
                new_ckt_segment = self.circuit

                self.final_mapping_ic(qiskit_ckt_object=new_ckt_segment)
                final_map = self.final_map
                self.set_initial_layout(final_map)

                new_ckt_segment.remove_final_measurements(inplace=True)
                IncrC_qc = IncrC_qc + new_ckt_segment

            beta = Parameter('b{}'.format(p))
            for node in range(logical_n):
                IncrC_qc.rx(beta,self.initial_layout[node])

            IncrC_qc.barrier()

        for i in range(logical_n):
            IncrC_qc.measure(self.initial_layout[i],i)

        self.circuit = IncrC_qc



    def instruction_parallelization(self,Sorted_Edges=None,Single_Layer=False):
        """
        This method is used for instruction parallelization.
        """
        Logical_qubits = self.zz_graph.nodes()
        if Sorted_Edges:
            Remaining_edges = Sorted_Edges.copy()
        else:
            Remaining_edges = list(self.zz_graph.edges())

        Current_layer = 'L0'
        Layer_occupancy = {Current_layer: list()}
        for qubit in Logical_qubits:
            Layer_occupancy[Current_layer].insert(len(Layer_occupancy[Current_layer]),[qubit,'FREE'])
        self.layer_zz_assignments[Current_layer] = list()

        while True:
            Unallocated_edges = list()
            allocated_op_count_in_this_layer = 0
            for edge in Remaining_edges:
                if allocated_op_count_in_this_layer >= self.Packing_Limit:
                    Unallocated_edges.insert(len(Unallocated_edges),edge)
                    continue

                n1, n2 = edge[0], edge[1]
                Free_among_the_two = 0

                for Occupancy_info_list in Layer_occupancy[Current_layer]:
                    if Occupancy_info_list[0] in edge:
                        if Occupancy_info_list[1] == 'OCCUPIED':
                            Unallocated_edges.insert(len(Unallocated_edges),edge)
                            break
                        else:
                            Free_among_the_two = Free_among_the_two + 1
                    if Free_among_the_two == 2:
                        n1_indx = Layer_occupancy[Current_layer].index([n1,'FREE'])
                        n2_indx = Layer_occupancy[Current_layer].index([n2,'FREE'])
                        Layer_occupancy[Current_layer][n1_indx] = [n1,'OCCUPIED']
                        Layer_occupancy[Current_layer][n2_indx] = [n2,'OCCUPIED']
                        self.layer_zz_assignments[Current_layer].insert(len(Layer_occupancy[Current_layer]),edge)
                        allocated_op_count_in_this_layer = allocated_op_count_in_this_layer + 1
                        break

            Remaining_edges = Unallocated_edges

            if Single_Layer:
                #print('Single layer formed!')
                self.layer_zz_assignments['R'] = list()
                for edge in Remaining_edges:
                    self.layer_zz_assignments['R'].insert(0,edge)
                break
            elif len(Remaining_edges) != 0:
                Next_layer = int(Current_layer[1:])+1
                Current_layer = 'L' + str(Next_layer)
                Layer_occupancy[Current_layer] = list()
                self.layer_zz_assignments[Current_layer] = list()
                for qubit in Logical_qubits:
                    Layer_occupancy[Current_layer].insert(len(Layer_occupancy[Current_layer]),[qubit,'FREE'])
            else:
                #print('All layers formed!')
                break

    def iterative_compilation(self):
        """
        This method is used for iterative compilation.
        """
        interchange = []
        LOO = []
        for L in self.layer_zz_assignments.keys():
            LOO.append(L)

        opt_target = 10e10
        opt_ckt = QuantumCircuit()

        while True:
            for layer_1 in range(len(LOO)):
                for layer_2 in range(layer_1+1,len(LOO)):
                    Temp = LOO.copy()
                    Temp[layer_1], Temp[layer_2] = Temp[layer_2], Temp[layer_1]

                    self.construct_circuit_iterc(LO=Temp)
                    trial_ckt = self.circuit

                    self.calc_cost(circ=trial_ckt,Target=self.IterC_Target)
                    trial_target = self.cost

                    if trial_target < opt_target:
                        interchange = [layer_1, layer_2]
                        opt_target = trial_target
                        opt_ckt = self.circuit


            if not interchange:
                self.circuit = opt_ckt
                break
            else:
                layer_1 = interchange[0]
                layer_2 = interchange[1]
                LOO[layer_1], LOO[layer_2] = LOO[layer_2], LOO[layer_1]
                #print('Interchanged: %s, %s, Cost: %s\n' % (layer_1, layer_2, opt_target))
                interchange = []

    def calc_cost(self,circ=None,Target='GC_2Q'):
        """
        This method is used to calculate cost of the compiled 
        circuit in terms of depth/2-qubit gate-count/estimated success probability.
        """
        if Target == 'GC_2Q':
            self.cost = circ.count_ops()[self.Native_2Q[0]]
        elif Target == 'D':
            self.cost = circ.depth()
        elif Target == 'ESP':
            self.circuit = circ
            self.estimate_sp()


    def estimate_sp(self):
        """
        This method estimates the success probability of a compiled circuit.
        """
        cir = self.circuit.copy()
        ESP = 1
        while True:
            if cir._data:
                k = cir._data.pop()
                GATE = k[0].__dict__['name']
                if GATE not in self.Basis_Gates:
                    continue
                QUB = []
                for i in range(len(k[1])):
                    QUB.append(k[1][i].index)
                if len(QUB) == 1:
                    ESP = ESP*float(self.QC_DATA[GATE][str(QUB[0])])
                else:
                    if '({},{})'.format(QUB[0],QUB[1]) in self.QC_DATA[GATE].keys():
                        ESP = ESP*float(self.QC_DATA[GATE]['({},{})'.format(QUB[0],QUB[1])])
                    elif '({},{})'.format(QUB[1],QUB[0]) in self.QC_DATA[GATE].keys():
                        ESP = ESP*float(self.QC_DATA[GATE]['({},{})'.format(QUB[1],QUB[0])])
                    else:
                        print('Please check the device configuration file for the following qubit-pair data: {}, {}'.format(QUB[0],QUB[1]))
            else:
                break
        self.cost = -math.log(ESP)

    def construct_circuit_iterc(self,LO=None):
        """
        This method constructs the circuit for iterative compilation.
        """
        n = len(self.zz_graph.nodes())
        qc = QuantumCircuit(n, n)

        # superposition state applying hadamard to all the qubits
        for node in self.zz_graph.nodes():
            qc.h(node)
        # change based on the mixing and phase separation layer architectures
        for p in range(1,self.Target_p+1):
            for L in LO:
                # phase seperation depends on the number of edges
                for edge in self.layer_zz_assignments[L]:
                    n1 = edge[0]
                    n2 = edge[1]
                    gamma = Parameter('g{}_{}_{}'.format(p,n1,n2))
                    qc.cx(n1, n2)
                    qc.rz(gamma, n2)
                    qc.cx(n1, n2)

            # mixing depends on the number of nodes rx gates
            beta = Parameter('b{}'.format(p))
            for node in self.zz_graph.nodes():
                qc.rx(beta,node)
            qc.barrier()
        qc.measure(range(n),range(n))
        trans_ckt = transpile(qc, coupling_map = self.Coupling_Map, \
                basis_gates = self.Basis_Gates, initial_layout = self.initial_layout, \
                optimization_level = self.Opt_Level, seed_transpiler = self.Trans_Seed, routing_method = self.Route_Method)
        self.circuit = trans_ckt


    def qasm_note(self,ckt=None):
        """
        This method prints notes on the compilation.
        """
        print('##################### Notes on the Output File #############################')
        if ckt:
            self.circuit = ckt
            self.estimate_sp()
            print('Depth: {}, Gate-count(2Q): {}, ESP: {}'.format(ckt.depth(), ckt.count_ops()[self.Native_2Q[0]], math.exp(-self.cost)))
            print('The circuit is written with beta/gamma parameters at different P-lavels (https://arxiv.org/pdf/1411.4028.pdf)')
            print('bX --> beta parameter at p=X')
            print('gX_Y_Z --> gamma parameter at p=X between *logical* qubit Y and Z. For the MaxCut problem of unweighted graphs, gX_Y1_Z1 = gX_Y2_Z2 (https://arxiv.org/pdf/1411.4028.pdf)')
        else: print('Compilation Error! Please check the input files.')
        print('############################################################################')


    def run_IP(self):
        """
        This method runs instruction parallelization.
        """
        self.instruction_parallelization()
        LO = self.layer_zz_assignments.keys()
        self.construct_circuit_iterc(LO)
        ckt = self.circuit
        ckt.qasm(filename=self.output_file_name)
        print('############################################################################')
        print('Instruction Parallelization-only Compilation (IP) completed! \nQASM File Written: {}'.format(self.output_file_name))
        self.qasm_note(ckt)


    def run_IterC(self,Target='D'):
        """
        This method runs iterative compilation.
        """
        self.set_iterc_target(Target)
        self.instruction_parallelization()
        self.iterative_compilation()
        ckt = self.circuit
        ckt.qasm(filename=self.output_file_name)
        print('############################################################################')
        print('Iterative Compilation (IterC) completed! \nQASM File Written: {}'.format(self.output_file_name))
        self.qasm_note(ckt=ckt)


    def run_IncrC(self,Variation_Aware=False):
        """
        This method runs incremental compilation.
        """
        self.set_incrc_var_awarenes(Variation_Aware)
        self.incremental_compilation()
        ckt = self.circuit
        ckt.qasm(filename=self.output_file_name)
        print('############################################################################')
        if Variation_Aware: print('Variation-aware Incremental Compilation (VIC) completed!\nQASM File Written: {}'.format(self.output_file_name))
        else: print('Incremental Compilation (IC) completed!\nQASM File Written: {}'.format(self.output_file_name))
        self.qasm_note(ckt)
