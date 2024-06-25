# -*- coding: utf-8 -*-
# Copyright 2018 Peter C Kroon
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Modifications by Christopher van den Berg on 21/06/2024
# Description of changes:
# * Refactored pysmiles.read_smiles.read_smiles() method logic into parse_smiles() method in order to use custom
# UndirectedMolecularGraph class
# * Directly used pysmiles.read_smiles._tokenize() private method logic to tokenize SMILES strings within
# parse_smiles()
# * Refactored pysmiles.smiles_helper.mark_aromatic_atoms() method logic into _set_aromatic_bond_orders() method to
# streamline for UndirectedMolecularGraph class use case

import logging
from queues import FifoQueue
from molecular_graphs import UndirectedMolecularGraph
from typing import List, Dict
from helpers.exception import InvalidSmilesError
from pysmiles import smiles_helper
from pysmiles.read_smiles import TokenType
from rdkit import Chem


class MolecularGraphSmilesParser:
    def smiles_to_molecular_graph(self, smiles_string: str):
        """
        Generates a molecular graph from a SMILES string.

        Parameters
        ----------
        smiles_string : str
            The SMILES string molecular representation of a given molecule.

        Returns
        -------
        UndirectedMolecularGraph
            An undirected molecular graph representation of the SMILES string.

        Raises
        ------
        InvalidSmilesError
            If the SMILES string is invalid.
        """
        if Chem.MolFromSmiles(smiles_string) is None:
            raise InvalidSmilesError()

        return self._parse_smiles(smiles_string)

    def _parse_smiles(self, smiles_string: str):
        """
       Parses a SMILES string to generate a molecular graph.

       This method was adapted from the original work by Peter C Kroon (pysmiles.read_smiles.read_smiles()).
       Changes include:
            * Refactored  method logic in order to use custom UndirectedMolecularGraph class
            * Attribute names made more verbose/descriptive

       Parameters
       ----------
       smiles_string : str
           The SMILES string representing the molecular structure.

       Returns
       -------
       UndirectedMolecularGraph
           An undirected molecular graph representation of the SMILES string.

       Raises
       ------
       ValueError
           If there is a conflict in bond orders or if a bond is specified incorrectly.
       KeyError
           If there are unmatched ring indices.
       InvalidSmilesError
           If an aromatic atom is detected outside a ring.
       """
        fifo_queue: FifoQueue = FifoQueue()
        molecular_graph: UndirectedMolecularGraph = UndirectedMolecularGraph(fifo_queue)
        bond_order: dict = {'-': 1, '=': 2, '#': 3, '$': 4, ':': 1.5, '.': 0}
        anchor = None
        node_id: int = 0
        default_bond_order: int = 1
        next_bond_order = None
        branches: List = []
        ring_nums: Dict = {}

        for token_type, token in self._tokenize(smiles_string):
            match token_type:
                case TokenType.ATOM:
                    molecular_graph.add_node(node_id, **smiles_helper.parse_atom(token))
                    if anchor is not None:
                        if next_bond_order is None:
                            next_bond_order = default_bond_order

                        molecular_graph.add_edge(anchor, node_id, bond_order=next_bond_order)
                        next_bond_order = None
                    anchor = node_id
                    node_id += 1
                case TokenType.BRANCH_START:
                    branches.append(anchor)
                case TokenType.BRANCH_END:
                    anchor = branches.pop()
                case TokenType.BOND_TYPE:
                    if next_bond_order is not None:
                        raise ValueError('Previous bond (order {}) not used. '
                                         'Overwritten by "{}"'.format(next_bond_order, token))
                    next_bond_order = bond_order[token]
                case TokenType.RING_NUM:
                    if token in ring_nums:
                        ring_node_id, ring_bond_order = ring_nums[token]
                        if next_bond_order is None and ring_bond_order is None:
                            next_bond_order = default_bond_order
                        elif bond_order is None:
                            next_bond_order = next_bond_order
                        elif next_bond_order is None:
                            next_bond_order = ring_bond_order
                        elif next_bond_order != ring_bond_order:
                            raise ValueError('Conflicting bond orders for ring '
                                             'between indices {}'.format(token))
                        if (node_id - 1, ring_node_id) in molecular_graph.edges:
                            raise ValueError('Edge specified by marker {} already '
                                             'exists'.format(token))
                        if node_id - 1 == ring_node_id:
                            raise ValueError('Marker {} specifies a bond between an '
                                             'atom and itself'.format(token))
                        if next_bond_order:
                            molecular_graph.add_edge(node_id - 1, ring_node_id, bond_order=next_bond_order)
                        next_bond_order = None
                        del ring_nums[token]
                    else:
                        if node_id == 0:
                            raise ValueError("Can't have a marker ({}) before an atom"
                                             "".format(token))
                        ring_nums[token] = (node_id - 1, next_bond_order)
                        next_bond_order = None
                case TokenType.EZSTEREO:
                    logging.getLogger(__name__).warning(
                        'E/Z stereochemical information, which is specified by "%s", will be discarded',
                        token
                    )

        if ring_nums:
            raise KeyError('Unmatched ring indices {}'.format(list(ring_nums.keys())))

        # Detect cycles in graph and assign aromatic bond order if any cycle atoms are aromatic
        cycles = molecular_graph.is_cyclic()

        if cycles:
            ring_node_ids = set()
            for cycle in cycles:
                ring_node_ids.update(cycle)
            non_ring_node_ids = set(molecular_graph.nodes) - ring_node_ids
            for node_id in non_ring_node_ids:
                if molecular_graph.nodes[node_id].get('aromatic', False):
                    raise InvalidSmilesError("An aromatic atom has been detected outside of a ring")

            self._set_aromatic_bond_orders(molecular_graph, cycles)

        return molecular_graph

    @staticmethod
    def _tokenize(smiles_str: str):
        """
        Iterates over a SMILES string and tokenizes it.

        This method was adapted from the original work by Peter C Kroon (pysmiles.read_smiles._tokenize()).
        Changes include:
            * No changes made (except for docstring)

        Parameters
        ----------
        smiles_str : str
            The SMILES string to iterate over.

        Yields
        ------
        tuple
            A tuple describing the type of token and the associated data.
        """
        organic_subset = 'B C N O P S F Cl Br I * b c n o s p'.split()
        smiles = iter(smiles_str)
        peek = None
        while True:
            char = peek if peek else next(smiles, '')
            peek = None
            if not char:
                break
            if char == '[':
                token = char
                for char in smiles:
                    token += char
                    if char == ']':
                        break
                yield TokenType.ATOM, token
            elif char in organic_subset:
                peek = next(smiles, '')
                if char + peek in organic_subset:
                    yield TokenType.ATOM, char + peek
                    peek = None
                else:
                    yield TokenType.ATOM, char
            elif char in '-=#$:.':
                yield TokenType.BOND_TYPE, char
            elif char == '(':
                yield TokenType.BRANCH_START, '('
            elif char == ')':
                yield TokenType.BRANCH_END, ')'
            elif char == '%':
                yield TokenType.RING_NUM, int(next(smiles, '') + next(smiles, ''))
            elif char in '/\\':
                yield TokenType.EZSTEREO, char
            elif char.isdigit():
                yield TokenType.RING_NUM, int(char)

    @staticmethod
    def _set_aromatic_bond_orders(molecular_graph: UndirectedMolecularGraph, cycles: List[List]):
        """
        Set all aromatic atom bond order to 1.5.

        This method was adapted from the original work by Peter C Kroon (pysmiles.smiles_helper.mark_aromatic_atoms()).
        Changes include:
            * Refactored method logic to streamline for UndirectedMolecularGraph class use case
            * Attribute names made more verbose/descriptive

        Parameters
        ----------
        molecular_graph : UndirectedMolecularGraph
            The undirected molecular graph.
        cycles : List[List]
            A list of cycles detected in the molecular graph.

        Returns
        -------
        None
        """
        for cycle in cycles:
            for i_node, j_node in molecular_graph.edges:
                if i_node not in cycle or j_node not in cycle:
                    continue
                if (molecular_graph.nodes[i_node].get('aromatic', False)
                        and molecular_graph.nodes[j_node].get('aromatic', False)):
                    molecular_graph.edges[i_node, j_node]['bond_order'] = 1.5
