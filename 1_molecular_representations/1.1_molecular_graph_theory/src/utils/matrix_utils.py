import pandas
import networkx


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
            A DataFrame where rows correspond to nodes and columns correspond to unique elements.
            Each cell contains 1 if the node has the corresponding element, otherwise 0 (i.e. will be one-hot encoded).

        Notes
        -----
        The DataFrame will have nodes as rows and unique attributes as columns, with a binary representation
        indicating the absence or presence of each element for each node.

        Examples
        --------
        >>> import networkx as nx
        >>> G = nx.Graph()
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
    # Extract all unique elements in the graph
    unique_elements = set(networkx.get_node_attributes(graph, node_attribute).values())
    unique_elements = sorted(unique_elements)  # Sort for consistent ordering

    # Initialize the matrix as a DataFrame
    num_nodes = len(graph.nodes)
    element_matrix = pandas.DataFrame(0, index=range(num_nodes), columns=unique_elements)

    # Populate the matrix
    for node, data in graph.nodes(data=True):
        element = data[node_attribute]
        element_matrix.at[node, element] = 1

    return element_matrix
