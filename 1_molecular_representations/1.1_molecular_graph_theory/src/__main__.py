#!/usr/bin/env python3
from queues import FifoQueue
from molecular_graphs import UndirectedMolecularGraph
from typing import Union, List
from smiles import MolecularGraphSmilesParser


def main():
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
    undirected_molecular_graph.add_edge("C", "D")
    undirected_molecular_graph.add_edge("G", "A")

    undirected_graph_a_c_path: Union[List, None] = undirected_molecular_graph.find_path("A", "C")
    undirected_graph_a_g_path: Union[List, None] = undirected_molecular_graph.find_path("A", "G")
    undirected_has_cycles: bool = undirected_molecular_graph.is_cyclic()


if __name__ == '__main__':
    main()
