"""
This module holds implementation of the initial qubit mapping strategies
for QAOA workloads described in the following articles:
    https://ieeexplore.ieee.org/document/9251960
    https://ieeexplore.ieee.org/document/9256490
"""
import networkx as nx

def create_initial_layout(coupling_graph, problem_zz_interactions_graph, \
    method = 'qaim', FW = None):
    """
    returns the initial layout as a list.
    Args:
        coupling_graph - hardware coupling graph
        problem_zz_interactions_graph - problem qaoa zz interactions graph
        method - chosen initial layout method (qaim/vqp)
        FW - qubit-to-qubit distances (floyd-warshall), if not provided,
            it is computed inside initial_layout method.
    Returns:
        initial layout in as a list.
    """
    layout = initial_layout(coupling_graph, problem_zz_interactions_graph, method = method, FW = FW)
    return list(layout.values())

def initial_layout(coupling_graph, problem_zz_interactions_graph, method = 'qaim', FW = None):
    """
    This method performs intelligent initial mapping (qaim or vqp) for qaoa workloads.
    """
    #calculate pairwise qubit-to-qubit distances, weights of the edges in the
    #coupling graph denotes 2-qubit gate success probabilities
    #to apply qaim, the weights should be 1 (qaim is variation-agnostic)
    if FW is None:
        FW = nx.floyd_warshall(coupling_graph)
    #initial empty layout assignment
    logical_to_physical_layout = {
        node: 'NA' for node in problem_zz_interactions_graph.nodes()
    }

    #sort program qubits based on usage
    sorted_program_qubits = sort_program_qubits(problem_zz_interactions_graph)

    #assign physical qubits to logical qubits using Noise-adaptive QAIM
    if method == 'vqp':
        qubit_strengths_dict = hardware_qubit_strength_vqp(coupling_graph)
    else:
        qubit_strengths_dict = hardware_qubit_strength_qaim(coupling_graph)

    sorted_physical_qubits = sort_physical_qubits(qubit_strengths_dict)
    allocated_physical_qubits = []

    for program_qubit in sorted_program_qubits:
        if placed_neighbors := [
            neigh
            for neigh in problem_zz_interactions_graph.neighbors(program_qubit)
            if logical_to_physical_layout[neigh] != 'NA'
        ]:
            unallocated_physical_neighbors = []
            for qubit in placed_neighbors:
                allocated_physical = logical_to_physical_layout[qubit]
                unallocated_physical_neighbors.extend(
                    element
                    for element in coupling_graph.neighbors(allocated_physical)
                    if element not in allocated_physical_qubits
                )

            #uniqify the unallocated physical neighbors
            unallocated_physical_neighbors = list(set(unallocated_physical_neighbors))

            #if none unallocated physical neighbors exists, pick the best available one
            if not unallocated_physical_neighbors:
                for physical_qubit in sorted_physical_qubits:
                    if physical_qubit not in allocated_physical_qubits:
                        logical_to_physical_layout[program_qubit] = physical_qubit
                        allocated_physical_qubits.append(physical_qubit)
                        break
                continue

            #sort the unallocated physical qubits based on
            #link strength / distance from placed neighbors - metric
            metric_of_unalloc_from_placed_neighbors = {}
            for element in unallocated_physical_neighbors:
                cumulative_distance = 0
                for qubit in placed_neighbors:
                    allocated_physical = logical_to_physical_layout[qubit]
                    cumulative_distance = cumulative_distance + FW[element][allocated_physical]
                metric_of_unalloc_from_placed_neighbors[element] = \
                    qubit_strengths_dict[element]/cumulative_distance
            #sorting procedure
            sorted_unallocated_physical_neighbors = []
            for element in unallocated_physical_neighbors:
                if not sorted_unallocated_physical_neighbors:
                    sorted_unallocated_physical_neighbors.append(element)
                    continue
                i = 0
                done = 'NO'
                for unc in sorted_unallocated_physical_neighbors:
                    if metric_of_unalloc_from_placed_neighbors[element] > \
                        metric_of_unalloc_from_placed_neighbors[unc]:
                        sorted_unallocated_physical_neighbors.insert(i,element)
                        done = 'YES'
                        break
                    i = i + 1
                if done == 'NO':
                    sorted_unallocated_physical_neighbors.append(element)

            #assign the best unallocated physical neighbor to the program qubit
            logical_to_physical_layout[program_qubit] = sorted_unallocated_physical_neighbors[0]
            allocated_physical_qubits.append(sorted_unallocated_physical_neighbors[0])

        else:
            for physical_qubit in sorted_physical_qubits:
                if physical_qubit not in allocated_physical_qubits:
                    logical_to_physical_layout[program_qubit] = physical_qubit
                    allocated_physical_qubits.append(physical_qubit)
                    break
    return logical_to_physical_layout

def sort_program_qubits(problem_zz_interactions_graph):
    """
    This method sorts the program qubits based on the
    number of zz gates per qubit.
    """
    problem_profile = problem_profiling(problem_zz_interactions_graph)
    sorted_program_qubits = []
    for key in problem_profile:
        if not sorted_program_qubits:
            sorted_program_qubits.append(key)
            #print(sorted_program_qubits)
            continue
        i = 0
        assigned = 'NO'
        for element in sorted_program_qubits:
            if problem_profile[element] < problem_profile[key]:
                sorted_program_qubits.insert(i,key)
                assigned = 'YES'
                break
            i = i + 1
        if assigned == 'NO':
            sorted_program_qubits.append(key)
            #print(sorted_program_qubits)
    return sorted_program_qubits

def problem_profiling(problem_zz_interactions_graph):
    """
    This method creates a program profile (zz gates/qubit).
    Args:
        ZZ interaction graph
    Returns:
        Dictionary (number of zz gates on each qubit)
    """
    pp_dict = {}
    for node in problem_zz_interactions_graph.nodes():
        try:
            pp_dict[node] = sum(1 for _ in problem_zz_interactions_graph.neighbors(node))
        except TypeError:
            pp_dict[node] = 0
            print(f'Unconnected node: {node}')
    return pp_dict

def hardware_qubit_strength_vqp(coupling_graph):
    """
    This method calculates qubit connectivity strengths from
    the harware coupling graph for vqp strategy.
    https://ieeexplore.ieee.org/document/9251960
    Args:
        hardware coupling graph (edge weights = 1)
    Returns:
         Dictionary holding the connectivity strenths of every qubit.
    """
    con_strengths = {}
    for node in coupling_graph.nodes():
        neighbours = coupling_graph.neighbors(node)
        strength = sum(coupling_graph[node][neigh]['weight'] for neigh in neighbours)
        con_strengths[node] = strength
    return con_strengths

def hardware_qubit_strength_qaim(coupling_graph):
    """
    This method calculates qubit connectivity strengths from
    the harware coupling graph for qaim strategy.
    https://ieeexplore.ieee.org/document/9256490
    Args:
        Hardware coupling graph (edge weights = success probability of native 2q gate)
    Returns:
        Dictionary holding the connectivity strenths of every qubit.
    """
    distances = nx.floyd_warshall(coupling_graph)
    con_strengths = {}
    for node in coupling_graph.nodes():
        con_strengths[node] = 0
        for node2 in coupling_graph.nodes():
            if distances[node][node2] in [1, 2]:
                con_strengths[node] += 1
    return con_strengths

def sort_physical_qubits(physical_qubit_profile):
    """
    This method sorts the physical qubits based on the their connectivity strengths
    https://ieeexplore.ieee.org/document/9251960
    https://ieeexplore.ieee.org/document/9256490
    Args:
        Dictionary holding the connectivity strenths of every qubit.
    Returns:
        Sorted list of the physical qubits.
    """
    sorted_physical_qubits = []
    for key in physical_qubit_profile.keys():
        if not sorted_physical_qubits:
            sorted_physical_qubits.append(key)
            #print(sorted_physical_qubits)
            continue
        i = 0
        assigned = 'NO'
        for element in sorted_physical_qubits:
            if physical_qubit_profile[element] < physical_qubit_profile[key]:
                sorted_physical_qubits.insert(i,key)
                assigned = 'YES'
                break
            i = i + 1
        if assigned == 'NO':
            sorted_physical_qubits.append(key)
            #print(sorted_physical_qubits)
    return sorted_physical_qubits

if __name__ == '__main__':
    #testing the above functions
    from qiskit.test.mock import FakeTokyo
    hw_coupling_graph = nx.Graph()
    device = FakeTokyo()
    cmap = device.configuration().coupling_map
    for init_layout_method in ['qaim', 'vqp']:
        for edge in cmap:
            if init_layout_method == 'vqp':
                weight = 1-device.properties().gate_error('cx', edge)
            else:
                weight = 1 #qaim
            hw_coupling_graph.add_edge(edge[0], edge[1], weight = weight)
        for seed in range(20):
            problem_zz_interaction_graph = nx.erdos_renyi_graph(19, 0.7, seed=seed)
            il = initial_layout(hw_coupling_graph, \
                problem_zz_interaction_graph, method = init_layout_method)
            assert len(set(il.keys())) == len(set(il.values()))
            print(il)
