#!/usr/bin/env python3
import networkx
import pandas
from typing import Tuple
from numpy import ndarray
from networkx import Graph
from smiles import MolecularGraphSmilesParser
from pysmiles.read_smiles import read_smiles
from rdkit.Chem import Mol, MolFromSmiles, GetAdjacencyMatrix
from utils import matrix_utils


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

    print(f'Ciprofloxacin Custom Molecular Graph = {custom_mol_graph.__str__()}\n')
    print(f'Node 0 Attribute Data = {custom_mol_graph.get_node_data(0)}\n')
    print(f'Edge (0, 1) Attribute Data = {custom_mol_graph.get_edge_data((0, 1))}\n')

    """ Part 2 """
    # Generate ciprofloxacin networkx.Graph object
    pysmiles_mol_graph: Graph = read_smiles(ciprofloxacin_smiles)
    print(f'Ciprofloxacin networkx.Graph from PySMILES = {pysmiles_mol_graph.__str__()}\n')

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

    print(f'{nx_adj_matrix}\n')

    # 2) Use RDKit to parse SMILES string into RDKit Mol object, convert to adjacency matrix and
    # obtain NetworkX Graph object
    ciprofloxacin_rdkit_mol = MolFromSmiles(ciprofloxacin_smiles)
    rdkit_adj_matrix: Mol = GetAdjacencyMatrix(ciprofloxacin_rdkit_mol, useBO=True)
    rdkit_mol_graph: Graph = networkx.from_numpy_array(rdkit_adj_matrix)

    print(f'{rdkit_adj_matrix}\n')

    """ Part 2 - Node Features """
    # Generate attribute matrix (tuple containing ndarray and matrix key) for `element` node attributes. Matrix shows
    # how many of each element-element bond types there are in the molecule (e.g. C-C, C-O, C-N etc.)
    pysmiles_element_matrix: Tuple = networkx.attr_matrix(
        pysmiles_mol_graph,
        node_attr='element'
    )
    print(f'{pysmiles_element_matrix}\n')

    # Generate one-hot encoded node attribute matrix for `element` node attribute
    one_hot_element_matrix = matrix_utils.node_attribute_matrix(pysmiles_mol_graph, 'element')
    print(f'{one_hot_element_matrix}\n')

    # Generate one-hot encoded node attribute matrix for `aromatic` node attribute
    one_hot_aromatic_matrix = matrix_utils.node_attribute_matrix(pysmiles_mol_graph, 'aromatic')
    print(f'{one_hot_aromatic_matrix}\n')

    # Generate one-hot encoded node attribute matrix for `charge` node attribute
    one_hot_charge_matrix = matrix_utils.node_attribute_matrix(pysmiles_mol_graph, 'charge')
    print(f'{one_hot_charge_matrix}\n')

    """ Part 3 - Edge Features Matrix """
    # Generate one-hot encoded edge attribute matrix for `bond_order` edge attribute
    one_hot_bond_order_matrix = matrix_utils.edge_attribute_matrix(pysmiles_mol_graph, 'order')
    print(f'{one_hot_bond_order_matrix}\n')

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
    print(f'Node 0 to Node 8 Path = {custom_mol_graph.find_path(0, 8)}\n')
    print(f'Ciprofloxacin Molecular Graph Connected Nodes = {custom_mol_graph.connected_components()}\n')

    """ Part 2 - Breadth-First Search Algorithm """
    print(f'Node 0 to Node 8 Shortest Path = {custom_mol_graph.find_shortest_path(0, 8)}\n')

    """ Part 3 - Traverse Graph and Detect Cycles """
    print(f'Ciprofloxacin Molecular Graph Cycles = {custom_mol_graph.is_cyclic()}\n')

    test = 'test'


if __name__ == '__main__':
    main()
