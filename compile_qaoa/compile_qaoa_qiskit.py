"""
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
"""
import math
import re
import os
from . initial_layout_qaoa import create_initial_layout
from random import shuffle
from qiskit.circuit import Parameter
from qiskit import QuantumCircuit, transpile, Aer, execute
import commentjson as json
import networkx as nx
from qiskit.quantum_info.analysis import hellinger_fidelity
from qiskit.converters import circuit_to_dag


class CompileQAOAQiskit:
    """
    This class implements the QAOA-specific compilation policies
    described in the following articles
    https://ieeexplore.ieee.org/document/9251960
    https://ieeexplore.ieee.org/document/9218558
    https://ieeexplore.ieee.org/document/9256490
    After crating the object, compilation can be performed with the chosen compilation
    policies using the following public methods:
        run_ip
        run_iter_c
        run_incr_c
    The current implementation only supports compilation with qiskit compiler backend.
    Necessary instrcutions are given under __compile_with_backend method docstring
    to add support for other compilers (e.g. tket).
    """
    def __init__(self, circuit_json = None,
    qc_json = None, config_json = None, out_circuit_file_name = 'QAOA.qasm'):
        """
        This method initializes necessary config variables.
        """
        self.supported_backends = ['qiskit']
        self.supported_initial_layout_strategies = ['qaim','vqp','random']
        self.__load_config(config_json)
        self.output_file_name = out_circuit_file_name
        self.__extract_qc_data(qc_json)
        self.layer_zz_assignments = {}
        self.zz_graph = self.qaoa_zz_graph(circuit_json)
        with open(circuit_json) as f:
            self.zz_dict = json.loads(f.read())
        self.initial_layout = list(range(len(self.zz_graph.nodes())))
        self.circuit = None
        self.sorted_ops = None
        self.cost = 10e10
        self.final_map = []
        [self.uncompiled_ckt, self.naive_ckt] = self.__naive_compilation()

    def __naive_compilation(self):
        """
        This method constructs performs a naive compilation with the qiskit backend.
        (gates are randomly ordered).
        """
        n = len(self.zz_graph.nodes())
        qc = QuantumCircuit(n, n)
        for node in self.zz_graph.nodes():
            qc.h(node)
        for p in range(1,self.Target_p+1):
            for zz in self.zz_graph.edges():
                n1 = zz[0]
                n2 = zz[1]
                gamma = Parameter('g{}_{}_{}'.format(p, n1, n2))
                qc.cx(n1, n2)
                qc.rz(gamma, n2)
                qc.cx(n1, n2)
            beta = Parameter('b{}'.format(p))
            for node in self.zz_graph.nodes():
                qc.rx(beta,node)
            qc.barrier()
        qc.measure(range(n), range(n))

        trans_ckt = self.__compile_with_backend(ckt_qiskit = qc)
        #qc = self.__fix_parameters(qc)
        filename = 'uncompiled_' + self.output_file_name
        qc.qasm(filename = filename)
        self.__fix_param_names(filename)
        filename = 'naive_compiled_' + self.output_file_name
        trans_ckt.qasm(filename = filename)
        self.__fix_param_names(filename)
        return [qc, trans_ckt]
    
    def __fix_parameters(self, circ):
        n = circ.num_qubits
        new_ckt = QuantumCircuit(n, n)
        dag = circuit_to_dag(circ)
        for node in dag.topological_op_nodes():
            if node._op.is_parameterized():
                param = node._op.params[0].name
                captures = re.search('(g\d+)_(\d+)_(\d+)', param).groups()
                if captures:
                    g = captures[0]
                    n1 = captures[1]
                    n2 = captures[2]
                    gamma = Parameter(g)
                    if '({}, {})'.format(n1, n2) in self.zz_dict.keys():
                        coeff = float(self.zz_dict['({}, {})'.format(n1, n2)])
                    else:
                        coeff = float(self.zz_dict['({}, {})'.format(n2, n1)])
                    new_op = node.op.copy()
                    new_op.params = [2*coeff*gamma]
                    try:
                        new_ckt.append(new_op, qargs = node.qargs)
                    except:
                        import pdb
                        pdb.set_trace()
                else:
                    beta = Parameter(param)
                    new_op = node.op.copy()
                    new_op.params = beta
                    new_ckt.append(new_op, qargs = node.qargs)
            else:
                new_ckt.append(node.op, qargs = node.qargs)
        return new_ckt
                
                
    def __load_config(self, config_json = None):
        """
        This method loads the variables in the config json file.
        """
        with open(config_json) as f:
            self.config = json.load(f)

        if 'Target_p' in self.config.keys():
            self.Target_p = int(self.config['Target_p'])
        else:
            self.Target_p = 1
        if 'Packing_Limit' in self.config.keys():
            self.Packing_Limit = float(self.config['Packing_Limit'])
        else:
            self.Packing_Limit = 10e10
        if 'Route_Method' in self.config.keys():
            self.Route_Method = self.config['Route_Method']
        else:
            self.Route_Method = 'sabre'
        if 'Trans_Seed' in self.config.keys():
            self.Trans_Seed = int(self.config['Trans_Seed'])
        else:
            self.Trans_Seed = 0
        if 'Opt_Level' in self.config.keys():
            self.Opt_Level = int(self.config['Opt_Level'])
        else:
            self.Opt_Level = 1
        if 'Backend' in self.config.keys():
            self.Backend = str(self.config['Backend'])
            assert self.Backend in self.supported_backends
        else:
            self.Backend = 'qiskit'
            assert self.Backend in self.supported_backends

    def __extract_qc_data(self, qc_file = None):
        """
        This method extracts hardware information from the QC json file.
        """
        with open(qc_file) as f:
            self.qc_data = json.load(f)

        try:
            self.native_2q = eval(self.qc_data['2Q'].strip('"'))
        except:
            self.native_2q = self.qc_data['2Q']
        try:
            self.native_1q = eval(self.qc_data['1Q'].strip('"'))
        except:
            self.native_1q = self.qc_data['1Q']

        self.basis_gates = self.native_2q + self.native_1q
        self.coupling_map = []
        for key in self.qc_data[str(self.native_2q[0])].keys():
            n1, n2 = eval(key)[0], eval(key)[1]
            if [n1, n2] not in self.coupling_map:
                self.coupling_map.append([n1, n2])
            if [n2, n1] not in self.coupling_map:
                self.coupling_map.append([n2, n1])
        self.__calc_qq_distances()

    def __calc_qq_distances(self):
        """
        This method calculates pairwise qubit-qubit distances using the floyd_warshall algorithm.
        """
        self.unweighted_undirected_coupling_graph = nx.Graph()
        self.weighted_undirected_coupling_graph = nx.Graph()
        for key, value in self.qc_data[str(self.native_2q[0])].items():
            n1, n2, sp = eval(key)[0], eval(key)[1], float(value)
            self.unweighted_undirected_coupling_graph.add_edge(n1, n2)
            self.weighted_undirected_coupling_graph.add_edge(n1, n2, weight=1/sp)
        self.qq_distances = nx.floyd_warshall(self.unweighted_undirected_coupling_graph)
        self.noise_aware_qq_distances = nx.floyd_warshall(self.weighted_undirected_coupling_graph)

    def __set_iter_c_target(self, target = 'GC_2Q'):
        """
        This method can be used to set the target of iterative compilation.
        """
        self.iter_c_target = target

    def __set_incrc_var_awareness(self, variation_aware = False):
        """
        This method can be used to set variation awareness in incremental compilation.
        """
        self.incr_c_var_awareness = variation_aware

    @staticmethod
    def qaoa_zz_graph(circ_json = None):
        """
        This method is used to create the MaxCut graph from the json file.
        """
        with open(circ_json) as f:
            data = json.loads(f.read())
        zz_graph = nx.Graph()
        for key, val in data.items():
            nodes = eval(key)
            zz_graph.add_edge(nodes[0], nodes[1])
        return zz_graph

    def __final_mapping_ic(self, qiskit_ckt_object):
        """
        This method finds the output logical-to-physical qubit mapping in a compiled circuit block.
        """
        qiskit_ckt_object.qasm(filename='tempfile')
        os.system('grep measure tempfile | awk \'{print $2, $4}\' > temp')
        qreg_creg_map = open('temp','r').readlines()

        fmap = {}
        for line in qreg_creg_map:
            elements = line.split(' ')
            physical_qs = elements[0]
            logical_qs = str(elements[1]).split(';')[0]
            physical_q = re.search('\[.*\]',physical_qs).group()
            physical_q = re.sub('\[|\]','',physical_q)
            logical_q = re.search('\[.*\]',logical_qs).group()
            logical_q = re.sub('\[|\]','',logical_q)
            fmap[logical_q[0:]] = int(physical_q[0:])

        final_map = []
        for i in range(qiskit_ckt_object.width() - qiskit_ckt_object.num_qubits):
            final_map.append(fmap[str(i)])
        os.system('rm temp')
        os.system('rm tempfile')
        self.final_map = final_map

    def __set_initial_layout(self, target_layout = None, initial_layout_method = None):
        """
        This method is used to set initial layout before starting compilation of any circuit block.
        """
        if target_layout:
            self.initial_layout = target_layout
        elif initial_layout_method:
            if initial_layout_method == 'qaim':
                self.initial_layout = create_initial_layout(self.weighted_undirected_coupling_graph, self.zz_graph, method = 'qaim')
            elif initial_layout_method == 'vqp':
                self.initial_layout = create_initial_layout(self.weighted_undirected_coupling_graph, self.zz_graph, method = 'vqp')
            elif initial_layout_method == 'random':
                self.initial_layout = list(range(len(self.zz_graph.nodes())))
                shuffle(self.initial_layout) #default: random
            else:
                raise ValueError

    def __sort_zz_by_qq_distances(self, unsorted_ops = None):
        """
        This method sorts the ZZ operations based on the control-target distances for the current mapping.
        """
        sorted_ops = []
        swap_distances_ops = {}

        for op in unsorted_ops:
            _physical_q1 = self.initial_layout[op[0]]
            _physical_q2 = self.initial_layout[op[1]]
            if not self.incr_c_var_awareness:
                swap_dist = self.qq_distances[_physical_q1][_physical_q2]
            else:
                swap_dist = self.noise_aware_qq_distances[_physical_q1][_physical_q2]
            swap_distances_ops[op] = swap_dist

        for op in unsorted_ops:
            if not sorted_ops:
                sorted_ops.append(op)
                continue
            i = 0
            for sop in sorted_ops:
                if swap_distances_ops[op] < swap_distances_ops[sop]:
                    sorted_ops.insert(i, op)
                    break
                i = i + 1
            if i == len(sorted_ops):
                sorted_ops.append(op)
        self.sorted_ops = sorted_ops

    def __construct_single_layer_ckt_ic(self, p):
        """
        This method constructs a single layer of the circuit in incremental compilation.
        """
        n = len(self.zz_graph.nodes())
        qc = QuantumCircuit(n, n)
        for zz in self.layer_zz_assignments['L0']:
            n1 = zz[0]
            n2 = zz[1]
            gamma = Parameter('g{}_{}_{}'.format(p, n1, n2))
            qc.cx(n1, n2)
            qc.rz(gamma, n2)
            qc.cx(n1, n2)
        qc.measure(range(n), range(n))
        trans_ckt = self.__compile_with_backend(ckt_qiskit = qc)
        self.circuit =  trans_ckt

    def __compile_with_backend(self, ckt_qiskit = None):
        """
        This method performs full/partial circuit compilation using the chosen backend.
        This method can be extended to support other compilers (e.g. tket).
        1) The target backend should support parametric circuit compilation from a given
        qiskit QuantumCircuit object with the target basis gates, given initial layout,
        and target hardware coupling map.
        initial layout format: [1 2 0 ...] (means q0 --> p1, q1 --> p2, q2 --> p0 ...)
                                qX: logical qubit, pX: physical qubit
        coupling map format: [[0,1],[1,0]...] (means the native 2q gate supported between
                                0 and 1, 1 and 0, ....)
        basis gates format: ['x', 'cx' ...]
        2) It should be able to recognize CX, RZ(theta), RX(theta), and H operations.
        3) It should return the compiled circuit as a qiskit QuantumCircuit object. If it
        does not do so, please convert the compiled circuit to a qiskit QuantumCircuit object.
        4) If you are adding a new backend, please update the supported_backends variable
        under __init__ as well.
        5) Use/add config variables in the Config.json (see under examples) file based on the
        supported features of the new backend.
        6) Update the __load_config method as well if you are adding new variables in Config.json.
        """
        assert isinstance(ckt_qiskit, QuantumCircuit)
        if self.Backend == 'qiskit':
            return transpile(ckt_qiskit, coupling_map = self.coupling_map,
                basis_gates = self.basis_gates, initial_layout = self.initial_layout,
                optimization_level = self.Opt_Level, seed_transpiler = self.Trans_Seed,
                routing_method = self.Route_Method)

    def __incremental_compilation(self):
        """
        This method is used for incremental compilation.
        """
        logical_n = len(self.zz_graph.nodes())
        physical_n = len(self.unweighted_undirected_coupling_graph.nodes())
        incr_c_qc = QuantumCircuit(physical_n, logical_n)

        for i in range(logical_n):
            incr_c_qc.h(self.initial_layout[i])

        for p in range(1,self.Target_p+1):
            remaining_ops = self.zz_graph.edges()
            while remaining_ops:
                self.__sort_zz_by_qq_distances(unsorted_ops = remaining_ops)
                sorted_ops = self.sorted_ops
                self.__instruction_parallelization(sorted_ops, single_layer = True)

                remaining_ops = self.layer_zz_assignments['R']

                self.__construct_single_layer_ckt_ic(p)
                new_ckt_segment = self.circuit

                self.__final_mapping_ic(qiskit_ckt_object = new_ckt_segment)
                final_map = self.final_map
                self.__set_initial_layout(final_map)

                new_ckt_segment.remove_final_measurements(inplace=True)
                incr_c_qc = incr_c_qc + new_ckt_segment

            beta = Parameter('b{}'.format(p))
            for node in range(logical_n):
                incr_c_qc.rx(beta,self.initial_layout[node])
            incr_c_qc.barrier()

        for i in range(logical_n):
            incr_c_qc.measure(self.initial_layout[i], i)

        self.circuit = incr_c_qc

    def __instruction_parallelization(self, sorted_edges = None, single_layer = False):
        """
        This method is used for instruction parallelization.
        """
        logical_qubits = self.zz_graph.nodes()
        if sorted_edges:
            remaining_edges = sorted_edges.copy()
        else:
            remaining_edges = list(self.zz_graph.edges())
        current_layer = 'L0'
        layer_occupancy = {current_layer: list()}
        for qubit in logical_qubits:
            layer_occupancy[current_layer].insert(len(layer_occupancy[current_layer]), [qubit,'FREE'])
        self.layer_zz_assignments[current_layer] = list()

        while True:
            unallocated_edges = list()
            allocated_op_count_in_this_layer = 0
            for edge in remaining_edges:
                if allocated_op_count_in_this_layer >= self.Packing_Limit:
                    unallocated_edges.insert(len(unallocated_edges), edge)
                    continue

                n1, n2 = edge[0], edge[1]
                free_among_the_two = 0

                for occupancy_info_list in layer_occupancy[current_layer]:
                    if occupancy_info_list[0] in edge:
                        if occupancy_info_list[1] == 'OCCUPIED':
                            unallocated_edges.insert(len(unallocated_edges), edge)
                            break
                        free_among_the_two = free_among_the_two + 1

                    if free_among_the_two == 2:
                        n1_indx = layer_occupancy[current_layer].index([n1, 'FREE'])
                        n2_indx = layer_occupancy[current_layer].index([n2, 'FREE'])
                        layer_occupancy[current_layer][n1_indx] = [n1, 'OCCUPIED']
                        layer_occupancy[current_layer][n2_indx] = [n2, 'OCCUPIED']
                        self.layer_zz_assignments[current_layer].insert(len(layer_occupancy[current_layer]), edge)
                        allocated_op_count_in_this_layer = allocated_op_count_in_this_layer + 1
                        break

            remaining_edges = unallocated_edges
            if single_layer:
                #print('Single layer formed!')
                self.layer_zz_assignments['R'] = list()
                for edge in remaining_edges:
                    self.layer_zz_assignments['R'].insert(0, edge)
                break
            elif len(remaining_edges) != 0:
                next_layer = int(current_layer[1:]) + 1
                current_layer = 'L' + str(next_layer)
                layer_occupancy[current_layer] = list()
                self.layer_zz_assignments[current_layer] = list()
                for qubit in logical_qubits:
                    layer_occupancy[current_layer].insert(len(layer_occupancy[current_layer]), [qubit,'FREE'])
            else:
                #print('All layers formed!')
                break

    def __iterative_compilation(self):
        """
        This method is used for iterative compilation.
        """
        interchange = []
        layer_order = []
        for l in self.layer_zz_assignments.keys():
            layer_order.append(l)

        opt_target = 10e10
        opt_ckt = QuantumCircuit()
        while True:
            for layer_1 in range(len(layer_order)):
                for layer_2 in range(layer_1+1, len(layer_order)):
                    temp = layer_order.copy()
                    temp[layer_1], temp[layer_2] = temp[layer_2], temp[layer_1]
                    self.__construct_circuit_iterc(layer_order = temp)
                    trial_ckt = self.circuit
                    self.__calc_cost(circ = trial_ckt, target = self.iter_c_target)
                    trial_target = self.cost

                    if trial_target < opt_target:
                        interchange = [layer_1, layer_2]
                        opt_target = trial_target
                        opt_ckt = self.circuit

            if not interchange:
                self.circuit = opt_ckt
                break

            layer_1 = interchange[0]
            layer_2 = interchange[1]
            layer_order[layer_1], layer_order[layer_2] = layer_order[layer_2], layer_order[layer_1]
            #print('Interchanged: %s, %s, Cost: %s\n' % (layer_1, layer_2, opt_target))
            interchange = []

    def __calc_cost(self, circ = None, target = 'GC_2Q'):
        """
        This method is used to calculate cost of the compiled
        circuit in terms of depth/2-qubit gate-count/estimated success probability.
        """
        if target == 'GC_2Q':
            self.cost = circ.count_ops()[self.native_2q[0]]
        elif target == 'D':
            self.cost = circ.depth()
        elif target == 'ESP':
            self.circuit = circ
            self.__estimate_sp()

    def __estimate_sp(self):
        """
        This method estimates the success probability of a compiled circuit.
        """
        cir = self.circuit.copy()
        ESP = 1
        while True:
            if cir._data:
                k = cir._data.pop()
                gate = k[0].__dict__['name']
                if gate not in self.basis_gates:
                    continue
                qub = []
                for i in range(len(k[1])):
                    qub.append(k[1][i].index)
                if len(qub) == 1:
                    ESP = ESP*float(self.qc_data[gate][str(qub[0])])
                else:
                    if '({},{})'.format(qub[0],qub[1]) in self.qc_data[gate].keys():
                        ESP = ESP*float(self.qc_data[gate]['({},{})'.format(qub[0], qub[1])])
                    elif '({},{})'.format(qub[1],qub[0]) in self.qc_data[gate].keys():
                        ESP = ESP*float(self.qc_data[gate]['({},{})'.format(qub[1], qub[0])])
                    else:
                        print('Please check the device configuration' +
                        'file for the following qubit-pair data: {}, {}'.format(qub[0], qub[1]))
            else:
                break
        self.cost = -math.log(ESP)

    def __construct_circuit_iterc(self, layer_order = None):
        """
        This method constructs the circuit for iterative compilation.
        """
        n = len(self.zz_graph.nodes())
        qc = QuantumCircuit(n, n)
        # superposition state applying hadamard to all the qubits
        for node in self.zz_graph.nodes():
            qc.h(node)
        # change based on the mixing and phase separation layer architectures
        for p in range(1, self.Target_p+1):
            for l in layer_order:
                # phase seperation depends on the number of edges
                for edge in self.layer_zz_assignments[l]:
                    n1 = edge[0]
                    n2 = edge[1]
                    gamma = Parameter('g{}_{}_{}'.format(p, n1, n2))
                    qc.cx(n1, n2)
                    qc.rz(gamma, n2)
                    qc.cx(n1, n2)
            # mixing depends on the number of nodes rx gates
            beta = Parameter('b{}'.format(p))
            for node in self.zz_graph.nodes():
                qc.rx(beta, node)
            qc.barrier()
        qc.measure(range(n), range(n))

        trans_ckt = self.__compile_with_backend(ckt_qiskit = qc)
        self.circuit = trans_ckt

    def __approximate_equivalence(self, ckt = None):
        """
        This method checks (approximate) equivalence of the compiled circuit with the original one
        by comparing the output measurements (at a fixed value of all the parameters).
        """
        bind_dic1 = {}
        bind_dic2 = {}
        val = 1
        for param in ckt.parameters:
            bind_dic1[param] = val
        for param in self.uncompiled_ckt.parameters:
            bind_dic2[param] = val

        ckt1 = ckt.bind_parameters(bind_dic1)
        ckt2 = self.uncompiled_ckt.bind_parameters(bind_dic2)
        backend_sim = Aer.get_backend('qasm_simulator')
        job_sim = execute([ckt1, ckt2], backend_sim, shots=1000000)
        result_sim = job_sim.result()
        counts1 = result_sim.get_counts(ckt1)
        counts2 = result_sim.get_counts(ckt2)
        return hellinger_fidelity(counts1, counts2) > 0.9

    def __qasm_note(self, ckt = None, pol = None):
        """
        This method prints notes on the compilation.
        """
        #optional, will not work for larger circuits due to finite sampling errors
        assert self.__approximate_equivalence(ckt)
        print('##################### Notes on the Output File #############################')

        if ckt:
            self.circuit = self.naive_ckt
            self.__estimate_sp()
            print('(naive) Depth: {}, gate-count(2Q): {}, ESP: {}'.format(self.naive_ckt.depth(),
            self.naive_ckt.count_ops()[self.native_2q[0]], math.exp(-self.cost)))

            self.circuit = ckt
            self.__estimate_sp()
            print('({}) Depth: {}, gate-count(2Q): {}, ESP: {}'.format(pol, ckt.depth(),
            ckt.count_ops()[self.native_2q[0]], math.exp(-self.cost)))

            print('The circuit is written with beta/gamma parameters' +
            'at different p-lavels (https://arxiv.org/pdf/1411.4028.pdf)')
            print('bX --> beta parameter at p=X')
            print('gX_Y_Z --> gamma parameter at p=X between *logical* ' +
            'qubit Y and Z. For the MaxCut problem of unweighted graphs, ' +
            'gX_Y1_Z1 = gX_Y2_Z2 (https://arxiv.org/pdf/1411.4028.pdf)')
        else:
            print('Compilation Error! Please check the input files.')

        print('############################################################################')

    def run_ip(self, initial_layout_method = 'qaim'):
        """
        This public method runs instruction parallelization and writes the output circuits
        in qasm format.
        args:
            No arguments required.
        """
        self.__set_initial_layout(initial_layout_method = initial_layout_method)
        self.__instruction_parallelization()
        layer_order = self.layer_zz_assignments.keys()
        self.__construct_circuit_iterc(layer_order)
        ckt = self.circuit
        filename = 'IP_' + self.output_file_name
        ckt.qasm(filename = filename)
        self.__fix_param_names(filename)

        print('############################################################################')
        print('Instruction Parallelization-only Compilation (IP) completed (initial layout: {})!'.format(initial_layout_method) +
        '\nQASM File Written: {}'.format('IP_' + self.output_file_name))
        self.__qasm_note(ckt, 'IP')

    def run_iter_c(self, target = 'D', initial_layout_method = 'qaim'):
        """
        This public method runs iterative compilation and writes the output circuits
        in qasm format.
        Args:
            Target minimization objective: D (depth), GC-2Q (two-qubit gate-count),
            ESP (estimated success probability)
        """
        self.__set_initial_layout(initial_layout_method = initial_layout_method)
        self.__set_iter_c_target(target)
        self.__instruction_parallelization()
        self.__iterative_compilation()
        ckt = self.circuit
        filename = 'IterC_' + self.output_file_name
        ckt.qasm(filename = filename)
        self.__fix_param_names(filename)
        print('############################################################################')
        print('Iterative Compilation (IterC) completed (initial layout: {})!'.format(initial_layout_method) +
        '\nQASM File Written: {}'.format('IterC_' + self.output_file_name))
        self.__qasm_note(ckt, 'IterC_' + target)

    def run_incr_c(self, variation_aware = False, initial_layout_method = 'qaim'):
        """
        This public method runs incremental compilation and writes the output circuits
        in qasm format.
        Args:
            variation_aware (boolean) - False to perform IC and True to perform VIC
        """
        self.__set_initial_layout(initial_layout_method = initial_layout_method)
        self.__set_incrc_var_awareness(variation_aware)
        self.__incremental_compilation()
        ckt = self.circuit

        print('############################################################################')
        if variation_aware:
            filename = 'VIC_' + self.output_file_name
            ckt.qasm(filename = filename)
            self.__fix_param_names(filename)
            print('Variation-aware Incremental Compilation (VIC) completed (initial layout: {})!'.format(initial_layout_method) +
            '\nQASM File Written: {}'.format('VIC_' + self.output_file_name))
            self.__qasm_note(ckt, 'VIC')
        else:
            filename = 'IC_' + self.output_file_name
            ckt.qasm(filename = filename)
            self.__fix_param_names(filename)
            print('Incremental Compilation (IC) completed (initial layout: {})!'.format(initial_layout_method) +
            '\nQASM File Written: {}'.format('IC_' + self.output_file_name))
            self.__qasm_note(ckt, 'IC')

    def __fix_param_names(self, filename):
        all_keys = self.zz_dict.keys()
        f = open(filename, 'r').readlines()
        out = open('{}_fixed'.format(filename), 'w')
        for line in f:
            captures = re.search('(g\d+)_(\d+)_(\d+)', line)
            if captures:
                captures = captures.groups()
                g = captures[0]
                n1 = captures[1]
                n2 = captures[2]
                if '({}, {})'.format(n1, n2) in all_keys:
                    coeff = 2*float(self.zz_dict['({}, {})'.format(n1, n2)])
                else:
                    coeff = 2*float(self.zz_dict['({}, {})'.format(n2, n1)])
                line = line.replace('{}_{}_{}'.format(g, n1, n2), str(coeff) + '*' + g)
            out.write(line)
        out.close()
        os.remove(filename)
        os.rename(filename + '_fixed', filename)



