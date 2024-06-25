import pandas
import networkx
import numpy
from typing import Set, List, Any


def node_attribute_matrix(graph: networkx.Graph, node_attribute: str):
    """
        Create a node attribute matrix that one-hot encodes a given node attribute.

        Parameters
        ----------
        graph : networkx.Graph
            The graph from which the node attributes are to be extracted.
        node_attribute : str
            The node attribute to be one-hot encoded and represented in the matrix.

        Returns
        -------
        pandas.DataFrame
            A DataFrame where rows correspond to nodes and columns correspond to a unique value of the specified node
            attribute. Each cell contains a binary value indicating absence (0) or presence (1) of the attribute value
            for each node (i.e. will be one-hot encoded).

        Examples
        --------
        >>> import networkx as nx
        >>> G = networkx.Graph()
        >>> G.add_node(0, element='C')
        >>> G.add_node(1, element='O')
        >>> G.add_node(2, element='C')
        >>> G.add_node(3, element='N')
        >>> G.add_edge(0, 1)
        >>> G.add_edge(1, 2)
        >>> G.add_edge(2, 3)
        >>> node_attribute_matrix(G, 'element')
           C  N  O
        0  1  0  0
        1  0  0  1
        2  1  0  0
        3  0  1  0
        """
    # Extract all unique `node_attribute` values in the graph
    unique_attr_values: Set = set(networkx.get_node_attributes(graph, node_attribute).values())

    # Sort for consistent ordering
    unique_attr_values: List = sorted(unique_attr_values)

    # Initialize the matrix as a DataFrame
    num_nodes: int = len(graph.nodes)
    node_attr_matrix: pandas.DataFrame = pandas.DataFrame(0, index=range(num_nodes), columns=unique_attr_values)

    # Populate the matrix
    for node, data in graph.nodes(data=True):
        attr_value = data[node_attribute]
        node_attr_matrix.at[node, attr_value] = 1

    return node_attr_matrix


def edge_attribute_matrix(graph: networkx.Graph, edge_attribute: str):
    """
    Create an edge attribute matrix that one-hot encodes a given edge attribute.

    Parameters
    ----------
    graph : networkx.Graph
        The input graph.
    edge_attribute : str
        The edge attribute to be represented in the matrix.

    Returns
    -------
    pandas.DataFrame
        A DataFrame where rows correspond to edges (represented as a tuple of nodes) and columns corresponds to a
        unique value of the specified edge attribute. Each cell contains a binary value indicating absence (0) or
        presence (1) of the attribute value for each edge (i.e. will be one-hot encoded).

    Examples
    --------
    >>> G = networkx.Graph()
    >>> G.add_edge(0, 1, bond_order=1)
    >>> G.add_edge(1, 2, bond_order=2)
    >>> G.add_edge(2, 3, bond_order=1)
    >>> edge_attribute_matrix(G, 'bond_order')
            1  2
    (0, 1)  1  0
    (1, 2)  0  1
    (2, 3)  1  0
    """
    unique_edge_attr_values: Set = set()

    # Extract all unique `edge_attribute` values in the graph
    for u_node, v_node, data in graph.edges(data=True):
        if edge_attribute in data:
            unique_edge_attr_values.add(data[edge_attribute])

    # Sort for consistent ordering
    unique_edge_attr_values: List = sorted(unique_edge_attr_values)

    # Initialize the matrix
    edge_keys: List = list(graph.edges)
    edge_attr_matrix: pandas.DataFrame = pandas.DataFrame(0, index=edge_keys, columns=unique_edge_attr_values)

    # Populate the matrix
    for u_node, v_node, data in graph.edges(data=True):
        if edge_attribute in data:
            attr_value: Any = data[edge_attribute]
            edge_attr_matrix.at[(u_node, v_node), attr_value] = 1

    return edge_attr_matrix
