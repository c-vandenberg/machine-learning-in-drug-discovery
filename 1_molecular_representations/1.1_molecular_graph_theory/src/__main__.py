#!/usr/bin/env python3
import os
import networkx
from typing import Union, Dict, Tuple, List, Set
from numpy import ndarray
from networkx import Graph
from smiles import MolecularGraphSmilesParser
from pysmiles.read_smiles import read_smiles
from rdkit.Chem import Mol, MolFromSmiles, GetAdjacencyMatrix
from utils import matrix_utils

os.environ['NUMEXPR_MAX_THREADS'] = '16'


def main():
    """
    1.1.1 Introduction To The Molecular Graph Representation

    Exercise:   Part 1 - Write custom MolecularGraphSmilesParser and UndirectedMolecularGraph classes to parse SMILES
                string of ciprofloxacin into an undirected molecular graph. Output node & edge attribute data

                Part 2 - In reality there is no need to create your own parser and molecular graph classes as many open
                source packages are available (e.g. `networkx`, `pysmiles`, `rdkit.Chem` etc.). Use PySMILES package to
                generate a NetworkX Graph object from ciprofloxacin SMILES string
    """
    """ Part 1 """
    # Instantiate SMILES to molecular graph parser
    molecular_graph_smiles_parser: MolecularGraphSmilesParser = MolecularGraphSmilesParser()
    ciprofloxacin_smiles: str = 'C1CNCCN1c(c2)c(F)cc3c2N(C4CC4)C=C(C3=O)C(=O)O'
    custom_mol_graph = molecular_graph_smiles_parser.smiles_to_molecular_graph(ciprofloxacin_smiles)

    # Output string representation of custom molecular graph
    custom_mol_graph_str_repr: str = custom_mol_graph.__str__()
    print(f'\033[1mCiprofloxacin Custom Molecular Graph:\033[0m\n{custom_mol_graph_str_repr}\n')

    # Output node 0 attribute data
    node_0_attr_data: Union[Dict, None] = custom_mol_graph.get_node_data(0)
    print(f'\033[1mNode 0 Attribute Data:\033[0m\n{node_0_attr_data}\n')

    # Output edge (0, 1) attribute data
    edge_0_1_attr_data: Union[Dict, None] = custom_mol_graph.get_edge_data((0, 1))
    print(f'\033[1mEdge (0, 1) Attribute Data:\033[0m\n{edge_0_1_attr_data}\n')

    """ Part 2 """
    # Generate ciprofloxacin networkx.Graph object and output string representation
    pysmiles_mol_graph: Graph = read_smiles(ciprofloxacin_smiles)
    pysmiles_mol_graph_str_repr: str = pysmiles_mol_graph.__str__()
    print(f'\033[1mCiprofloxacin networkx.Graph from PySMILES:\033[0m\n{pysmiles_mol_graph_str_repr}\n')

    """
    1.1.2 Mathematical Definition of a Graph
    
    Exercise:   Part 1 - Adjacency Matrix 
                    Use two approaches to obtain adjacency matrix:
                        1) Generate adjacency matrix from PySMILES NetworkX Graph object
                        2) Use RDKit to parse SMILES string into RDKit Mol object, convert to adjacency matrix and 
                        obtain NetworkX Graph object
                    
                Part 2 - Node Features Matrix
                    Use PySMILES molecular graph to construct node features matrix that utilises one-hot encoding
                    
                Part 3 - Edge Features Matrix
                    Use PySMILES molecular graph to construct edge features matrix that utilises one-hot encoding
    """

    """  Part 1 - Adjacency Matrix """
    # 1) Generate adjacency matrix from PySMILES NetworkX Graph object
    nx_adj_matrix: ndarray = networkx.adjacency_matrix(
        pysmiles_mol_graph
    ).todense()

    print(f'\033[1mCiprofloxacin PySMILES adjacency matrix:\033[0m\n{nx_adj_matrix}\n')

    # 2) Use RDKit to parse SMILES string into RDKit Mol object, convert to adjacency matrix and
    # obtain NetworkX Graph object
    ciprofloxacin_rdkit_mol = MolFromSmiles(ciprofloxacin_smiles)
    rdkit_adj_matrix: Mol = GetAdjacencyMatrix(ciprofloxacin_rdkit_mol, useBO=True)
    rdkit_mol_graph: Graph = networkx.from_numpy_array(rdkit_adj_matrix)

    print(f'\033[1mCiprofloxacin RDKit adjacency matrix:\033[0m\n{rdkit_adj_matrix}\n')

    """ Part 2 - Node Features """
    # Generate attribute matrix (tuple containing ndarray and matrix key) for `element` node attributes. Matrix shows
    # how many of each element-element bond types there are in the molecule (e.g. C-C, C-O, C-N etc.)
    pysmiles_element_matrix: Tuple = networkx.attr_matrix(
        pysmiles_mol_graph,
        node_attr='element'
    )
    print(f'\033[1mCiprofloxacin element bond types node attribute matrix:\033[0m\n{pysmiles_element_matrix}\n')

    # Generate one-hot encoded node attribute matrix for `element` node attribute
    one_hot_element_matrix = matrix_utils.node_attribute_matrix(pysmiles_mol_graph, 'element')
    print(f'\033[1mCiprofloxacin one-hot encoded `element` node attribute matrix:\033[0m\n{one_hot_element_matrix}\n')

    # Generate one-hot encoded node attribute matrix for `aromatic` node attribute
    one_hot_aromatic_matrix = matrix_utils.node_attribute_matrix(pysmiles_mol_graph, 'aromatic')
    print(f'\033[1mCiprofloxacin one-hot encoded `aromatic` node attribute matrix:\033[0m\n{one_hot_aromatic_matrix}\n')

    # Generate one-hot encoded node attribute matrix for `charge` node attribute
    one_hot_charge_matrix = matrix_utils.node_attribute_matrix(pysmiles_mol_graph, 'charge')
    print(f'\033[1mCiprofloxacin one-hot encoded `charge` node attribute matrix:\033[0m\n{one_hot_charge_matrix}\n')

    """ Part 3 - Edge Features Matrix """
    # Generate one-hot encoded edge attribute matrix for `order` (bond order) edge attribute
    one_hot_bond_order_matrix = matrix_utils.edge_attribute_matrix(pysmiles_mol_graph, 'order')
    print(f'\033[1mCiprofloxacin one-hot encoded `order` (bond order) edge attribute matrix:\033[0m\n'
          f'{one_hot_bond_order_matrix}\n')

    """
    1.1.3 Graph Traversal Algorithms

    Exercise:   Part 1 - Depth-First Search Algorithm
                    Add method to custom UndirectedMolecularGraph class that utilises a depth-first search (DFS) 
                    algorithm to traverse the graph. Utilise this DFS method to find a path between two nodes, and to 
                    detect all connected nodes in the graph.
                
                Part 2 - Breadth-First Search Algorithm
                    Add method to custom UndirectedMolecularGraph class that utilises a breadth-first search (BFS) 
                    algorithm to traverse the graph. Utilise this BFS method to find the shortest path between two 
                    nodes.
                
                Part 3 - Traverse Graph and Detect Cycles
                    Add method(s) to detect and return any cycles in the molecular graph.
    """

    """ Part 1 - Depth-First Search Algorithm """
    # Traverse graph using DFS to find path between two nodes
    node_0_node_8_path: Union[List, None] = custom_mol_graph.find_path(0, 8)
    print(f'\033[1mNode 0 to Node 8 Path:\033[0m\n{node_0_node_8_path}\n')

    # Traverse graph using DFS to output all connected nodes
    connected_nodes: Set = custom_mol_graph.connected_components()
    print(f'\033[1mCiprofloxacin Molecular Graph Connected Nodes:\033[0m\n{connected_nodes}\n')

    """ Part 2 - Breadth-First Search Algorithm """
    node_0_node_8_shortest_path: Union[List, None] = custom_mol_graph.find_shortest_path(0, 8)
    print(f'\033[1mNode 0 to Node 8 Shortest Path:\033[0m\n{node_0_node_8_shortest_path}\n')

    """ Part 3 - Traverse Graph and Detect Cycles """
    cycles: Union[List[List], None] = custom_mol_graph.is_cyclic()
    print(f'\033[1mCiprofloxacin Molecular Graph Cycles:\033[0m\n{cycles}\n')


if __name__ == '__main__':
    main()
