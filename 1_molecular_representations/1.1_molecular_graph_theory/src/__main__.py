#!/usr/bin/env python3
from queues import FifoQueue
from molecular_graphs import UndirectedMolecularGraph, DirectedMolecularGraph
from typing import Union, List, Set, Dict
from smiles import MolecularGraphSmilesParser


def main():
    # Instantiate SMILES to molecular graph parser
    molecular_graph_smiles_parser: MolecularGraphSmilesParser = MolecularGraphSmilesParser()

    ciprofloxacin_smiles: str = 'C1CNCCN1c(c2)c(F)cc3c2N(C4CC4)C=C(C3=O)C(=O)O'
    molecular_graph = molecular_graph_smiles_parser.smiles_to_molecular_graph(ciprofloxacin_smiles)
    node_0_data: Union[Dict, None] = molecular_graph.get_node_data(0)
    edge_0_1_data: Union[Dict, None] = molecular_graph.get_edge_data((0, 1))
    node_0_8_path = molecular_graph.find_path(0, 8)
    node_0_8_shortest_path = molecular_graph.find_shortest_path(0, 8)
    connected_nodes: Set = molecular_graph.connected_components()
    is_cyclic: Union[List, None] = molecular_graph.is_cyclic()

    print(molecular_graph.__str__())


if __name__ == '__main__':
    main()
