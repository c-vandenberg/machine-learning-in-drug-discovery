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
    Part 1: Use custom MolecularGraphSmilesParser and UndirectedMolecularGraph classes to parse SMILES string of
            ciprofloxacin into an undirected molecular graph. Test class method logic by printing returns.

            In reality there is no need to create your own parser and molecular graph classes as many open source
            packages available (e.g. `networkx`, `pysmiles`, `rdkit.Chem` etc.).
    """
    # Instantiate SMILES to molecular graph parser
    molecular_graph_smiles_parser: MolecularGraphSmilesParser = MolecularGraphSmilesParser()

    ciprofloxacin_smiles: str = 'C1CNCCN1c(c2)c(F)cc3c2N(C4CC4)C=C(C3=O)C(=O)O'
    custom_mol_graph = molecular_graph_smiles_parser.smiles_to_molecular_graph(ciprofloxacin_smiles)

    print(f'Ciprofloxacin Molecular = {custom_mol_graph.__str__()}\n')

    print(f'Node 0 Attribute Data = {custom_mol_graph.get_node_data(0)}\n')
    print(f'Edge (0, 1) Attribute Data = {custom_mol_graph.get_edge_data((0, 1))}\n')
    print(f'Node 0 to Node 8 Path = {custom_mol_graph.find_path(0, 8)}\n')
    print(f'Node 0 to Node 8 Shortest Path = {custom_mol_graph.find_shortest_path(0, 8)}\n')

    print(f'Ciprofloxacin Molecular Graph Connected Nodes = {custom_mol_graph.connected_components()}\n')
    print(f'Ciprofloxacin Molecular Graph Cycles = {custom_mol_graph.is_cyclic()}\n')

    """
    Part 2.1: Use PySMILES to parse SMILES string into NetworkX Graph object and convert to adjacency matrix
    """
    # Generate ciprofloxacin networkx.Graph object
    pysmiles_mol_graph: Graph = read_smiles(ciprofloxacin_smiles)

    # Print ciprofloxacin adjacency matrix
    nx_adj_matrix: ndarray = networkx.adjacency_matrix(
        pysmiles_mol_graph
    ).todense()

    print(f'{nx_adj_matrix}\n')

    """
    Part 2.2: Use RDKit to parse SMILES string into RDKit Mol object, convert to adjacency matrix and obtain NetworkX 
            Graph object
    """
    # Generate ciprofloxacin Mol object
    ciprofloxacin_rdkit_mol = MolFromSmiles(ciprofloxacin_smiles)

    # Generate ciprofloxacin adjacency matrix
    rdkit_adj_matrix: Mol = GetAdjacencyMatrix(ciprofloxacin_rdkit_mol, useBO=True)

    # Convert RDKit adjacency matrix to NetworkX Graph
    rdkit_mol_graph: Graph = networkx.from_numpy_array(rdkit_adj_matrix)

    print(f'{rdkit_mol_graph}\n')

    """
    Part 3: Use PySMILES molecular graph to construct node features matrices
    """
    # Generate attribute matrix (tuple containing ndarray and matrix key) for `element` node attributes. Matrix shows
    # how many of each element-element bond types there are in the molecule (e.g. C-C, C-O, C-N etc.)
    pysmiles_element_matrix: Tuple = networkx.attr_matrix(
        pysmiles_mol_graph,
        node_attr='element'
    )
    print(f'{pysmiles_element_matrix}\n')

    # Generate one-hot encoded node attribute matrix for `element` node attributes
    one_hot_element_matrix = matrix_utils.node_attribute_matrix(pysmiles_mol_graph, 'element')

    test = 'test'


if __name__ == '__main__':
    main()
