#!/usr/bin/env python3
from queues import FifoQueue
from molecular_graphs import UndirectedMolecularGraph, DirectedMolecularGraph
from typing import Union, List, Set
from smiles import MolecularGraphSmilesParser


def main():
    directed_fifo_queue: FifoQueue = FifoQueue()
    directed_molecular_graph: DirectedMolecularGraph = DirectedMolecularGraph(directed_fifo_queue)

    # Add nodes
    directed_molecular_graph.add_node("A")
    directed_molecular_graph.add_node("B")
    directed_molecular_graph.add_node("C")
    directed_molecular_graph.add_node("D")
    directed_molecular_graph.add_node("E")
    directed_molecular_graph.add_node("F")
    directed_molecular_graph.add_node("G")
    directed_molecular_graph.add_node("H")
    directed_molecular_graph.add_node("I")

    directed_molecular_graph.add_edge("A", "B")
    directed_molecular_graph.add_edge("B", "C")
    directed_molecular_graph.add_edge("B", "G")
    directed_molecular_graph.add_edge("C", "G")
    directed_molecular_graph.add_edge("G", "A")

    directed_graph_connected_nodes: Set = directed_molecular_graph.connected_components()
    directed_has_cycles: bool = directed_molecular_graph.is_cyclic()

    directed_graph_shortest_a_g_path: Union[List, None] = directed_molecular_graph.find_shortest_path(
        "A", "C"
    )

    # Create undirected graph
    undirected_fifo_queue: FifoQueue = FifoQueue()
    undirected_molecular_graph: UndirectedMolecularGraph = UndirectedMolecularGraph(undirected_fifo_queue)
    molecular_graph_smiles_parser: MolecularGraphSmilesParser = MolecularGraphSmilesParser()

    acetaminophen_smiles: str = 'CC(=O)Nc1ccc(cc1)O'
    molecular_graph = molecular_graph_smiles_parser.smiles_to_molecular_graph(acetaminophen_smiles)

    print(molecular_graph.__str__())

    # Add nodes
    undirected_molecular_graph.add_node("A")
    undirected_molecular_graph.add_node("B")
    undirected_molecular_graph.add_node("C")
    undirected_molecular_graph.add_node("D")
    undirected_molecular_graph.add_node("E")
    undirected_molecular_graph.add_node("F")
    undirected_molecular_graph.add_node("G")
    undirected_molecular_graph.add_node("H")
    undirected_molecular_graph.add_node("I")

    # Add edges
    undirected_molecular_graph.add_edge("A", "B")
    undirected_molecular_graph.add_edge("B", "C")
    undirected_molecular_graph.add_edge("B", "G")
    undirected_molecular_graph.add_edge("C", "G")
    undirected_molecular_graph.add_edge("G", "A")

    undirected_graph_a_c_path: Union[List, None] = undirected_molecular_graph.find_path("A", "C")
    undirected_graph_a_g_path: Union[List, None] = undirected_molecular_graph.find_path("A", "G")
    cycles = undirected_molecular_graph.is_cyclic()
    undirected_graph_a_g_shortest_path: Union[List, None] = undirected_molecular_graph.find_shortest_path("A", "G")
    connected_nodes = undirected_molecular_graph.connected_components()
    test = 'test'


if __name__ == '__main__':
    main()
