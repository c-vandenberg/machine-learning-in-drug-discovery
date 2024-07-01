import logging
from typing import Any, Union, List, Set, Dict, Tuple
from queues import FifoQueue
from helpers.exception import NonexistentNodeError, NonexistentEdgeError, CycleError

logger = logging.getLogger(__name__)


class BaseMolecularGraph:
    """
    Base class for a molecular graph data structure.

    Methods
    -------
    nodes() -> Dict[Union[int, float], List[Dict]]
        Getter for the graph nodes and their edges.
    add_node(node: Any)
        Add node to graph if it doesn't already exist.
    delete_node(node_to_delete: Any, default_return=None) -> None
        Delete a node from the graph if it exists, along with any edges that contain the node.
    get_node_data(node: Any, default_return=None) -> Any
        Get the attribute dictionary associated with a node if it exists.
    add_edge(from_node: Any, to_node: Any, weight: Union[int, float, None] = None)
        Add edge between two nodes in graph.
    delete_edge(edge_to_delete: Tuple[Any, Any]) -> None
        Delete an edge from the graph if it exists.
    get_edge_data(edge: Tuple, default_return=None) -> Any
        Get the attribute dictionary associated with an edge if it exists.
    find_path(start_node: Any, end_node: Union[None, Any] = None) -> Union[List, None]
        Find a path between two nodes using DFS (will not necessarily be the shortest path).
    find_shortest_path(start_node: Any, end_node: Union[None, Any] = None) -> Union[List, None]
        Find the shortest path between two nodes using BFS.
    connected_components() -> Set
        Return all connected components in the graph.
    """

    def __init__(self, fifo_queue: FifoQueue):
        self._fifo_queue = fifo_queue
        self._nodes = dict()
        self._edges = dict()

    @property
    def nodes(self) -> Dict:
        """
        Getter for the graph nodes

        Returns
        -------
        Dict
            A dictionary where keys are nodes and the values are sub dictionaries of the node attributes
        """
        return self._nodes

    @property
    def edges(self) -> Dict:
        """
        Getter for the graph edges

        Returns
        -------
        Dict
            A dictionary where keys are tuples representing the graph edges, and the values are sub dictionaries of the
            edge attributes
        """
        return self._edges

    def add_node(self, node: Any, **attr) -> None:
        """
        Add node to graph if it doesn't already exist.

        Parameters
        ----------
        node : Any
            The node to be added to the graph.
        **attr : Dict
            Any additional node attributes.

        Returns
        -------
        None
        """
        if node not in self._nodes:
            # Nodes are defined as a dictionary with 'neighbors' key and attributes
            self._nodes[node] = {'neighbours': [], **attr}
        else:
            # If node already exists, update the attributes
            self._nodes[node].update(attr)

    def delete_node(self, node_to_delete: any, default_return=None):
        """
        Delete a node from the graph if it exists, along with any edges that contain the node.

        Parameters
        ----------
        node_to_delete : Any
            The node to be deleted from the graph.
        default_return : Any, optional
            The value to return if the node does not exist. Defaults to None.

        Returns
        -------
        None

        Raises
        ------
        NonexistentNodeError
            If the node is not present in the graph.
        """
        try:
            self._validate_nodes(node_to_delete)

            # Remove node from neighbours list of all connected nodes
            for node in self._nodes:
                if node_to_delete in self._nodes[node]['neighbours']:
                    self._nodes[node]['neighbours'].remove(node_to_delete)

            # Remove node
            del self._nodes[node_to_delete]

            edges_to_remove: List = [edge for edge in self._edges if node_to_delete in edge]

            # Remove any edges that contain node
            for edge in edges_to_remove:
                del self._edges[edge]
        except NonexistentNodeError as e:
            logger.exception(e.message)

            return default_return

    def get_node_data(self, node: any, default_return=None) -> Union[Dict, None]:
        """
        Get the attribute dictionary associated with a node if it exists.

        Parameters
        ----------
        node : Any
            The node whose data is to be retrieved.
        default_return : Any, optional
            The value to return if the node does not exist. Defaults to None.

        Returns
        -------
        Union[Dict, None]
            The attribute dictionary of the node if it exists, else returns `default_return`.

        Raises
        ------
        NonexistentNodeError
            If the node is not present in the graph.
        """
        try:
            self._validate_nodes(node)

            return self._nodes[node]
        except NonexistentNodeError as e:
            logger.exception(e.message)

            return default_return

    def add_edge(self, from_node: Any, to_node: Any, **attr) -> None:
        """
        Add edge between two nodes in graph.

        Parameters
        ----------
        from_node : Any
            The starting node of the edge.
        to_node : Any
            The ending node of the edge.
        **attr : Dict
            Any additional edge attributes.

        Returns
        -------
        None
        """
        self._validate_nodes(from_node, to_node)
        self._add_edge(from_node, to_node, **attr)

    def delete_edge(self, edge_to_delete: Tuple) -> None:
        """
        Delete an edge from the graph if it exists.

        Parameters
        ----------
        edge_to_delete : Tuple[Any, Any]
            The edge to be deleted from the graph, represented as a tuple of two nodes.

        Returns
        -------
        None

        Raises
        ------
        NonexistentEdgeError
            If the edge is not present in the graph.
        """
        try:
            self._validate_edge(edge_to_delete)

            # Remove two edge nodes from each other's neighbours list if they are present (accounts for directed and
            # non-directed graphs)
            node_u, node_v = edge_to_delete
            if node_v in self._nodes[node_u]['neighbours']:
                self._nodes[node_u]['neighbours'].remove(node_v)
            if node_u in self._nodes[node_v]['neighbours']:
                self._nodes[node_v]['neighbours'].remove(node_u)

            # Remove edge
            del self._edges[edge_to_delete]
        except NonexistentEdgeError as e:
            logger.exception(e.message)

    def get_edge_data(self, edge: Tuple, default_return=None) -> Union[Dict, None]:
        """
        Get the attribute dictionary associated with an edge if it exists.

        Parameters
        ----------
        edge : Tuple
            The edge whose data is to be retrieved.
        default_return : Any, optional
            The value to return if the node does not exist. Defaults to None.

        Returns
        -------
        Union[Dict, None]
            The attribute dictionary of the edge if it exists, else returns `default_return`.

        Raises
        ------
        NonexistentEdgeError
            If the edge is not present in the graph.
        """
        try:
            self._validate_edge(edge)

            return self._edges[edge]
        except NonexistentEdgeError as e:
            logger.exception(e.message)

            return default_return

    def find_path(self, start_node: Any, end_node: Union[None, Any] = None) -> Union[List, None]:
        """
        Find a path between two nodes using DFS (will not necessarily be the shortest path).

        Parameters
        ----------
        start_node : Any
            The starting node.
        end_node : Union[None, Any]. Optional
            The target node, by default None.

        Returns
        -------
        Union[List, None]
            A list representing the path if it exists, else None.

        Raises
        ------
        NonexistentNodeError
            If either the start node or the end node is not present in the graph.
        """
        try:
            self._validate_nodes(start_node, end_node)

            predecessor = {start_node: None}
            visited_nodes = set()
            path = self._dfs(start_node, end_node, predecessor, visited_nodes)

            return self._reconstruct_path(predecessor, end_node) if path else None
        except NonexistentNodeError as e:
            logger.exception(e.message)

            return None

    def find_shortest_path(self, start_node: Any, end_node: Union[None, Any] = None) -> Union[List, None]:
        """
        Find the shortest path between two nodes using BFS.

        Parameters
        ----------
        start_node : Any
            The starting node.
        end_node : Union[None, Any]. Optional
            The target node, by default None.

        Returns
        -------
        Union[List, None]
            A list representing the shortest path if it exists, else None.

        Raises
        ------
        NonexistentNodeError
            If either the start node or the end node is not present in the graph.
        """
        try:
            self._validate_nodes(start_node, end_node)

            return self._bfs(start_node, end_node)
        except NonexistentNodeError as e:
            logger.exception(e.message)

            return None

    def connected_components(self) -> Set:
        """
        Return all connected components in the graph using depth-first search (DFS) algorithm to traverse graph.

        Returns
        -------
        set
            A set of nodes that have neighbors (connected components).
        """
        visited_nodes = set()

        for node in self._nodes:
            if self._has_traversable_neighbours(self._nodes[node]) and node not in visited_nodes:
                self._dfs_connected_component(node, visited_nodes)

        return visited_nodes

    def __str__(self) -> str:
        """
        Return the string representation of the graph.

        Returns
        -------
        str
            A string representation of the nodes and their edges.
        """
        return str(self._nodes)

    def _add_edge(self, from_node: Any, to_node: Any, **attr) -> None:
        """
        Add edge between two nodes in graph (to be implemented/overridden by subclasses).

        Parameters
        ----------
        from_node : Any
            The starting node of the edge.
        to_node : Any
            The ending node of the edge.
        **attr : Dict
            Any additional edge attributes.

        Returns
        -------
        None
        """
        pass

    def _validate_nodes(self, start_node: Any, end_node: Union[None, Any] = None) -> None:
        """
        Validate that the nodes exist in the graph.

        Parameters
        ----------
        start_node : Any
            The starting node.
        end_node : Union[None, Any]. Optional
            The ending node, by default None.

        Returns
        -------
        None

        Raises
        ------
        NonexistentNodeError
            If either the start node or the end node is not present in the graph.
        """
        if start_node not in self._nodes or (end_node not in self._nodes and end_node):
            missing_node: Any = start_node if start_node not in self._nodes else end_node

            raise NonexistentNodeError(f'Node {missing_node} not present in Graph')

    def _validate_edge(self, edge: Tuple) -> None:
        """
        Validate that the edge exist in the graph.

        Parameters
        ----------
        edge : Tuple
            The edge to validate.

        Returns
        -------
        None

        Raises
        ------
        NonexistentEdgeError
            If either the edge is not present in the graph.
        """
        if edge not in self._edges:
            edge_str: str = ' -> '.join(map(str, edge))
            raise NonexistentEdgeError(f'Edge {edge_str} not present in Graph')

    def _detect_cycles(self, current_node: Any, visited: Dict[Any, bool], path: List, cycles: List,
                       parent_node: Any = None, is_undirected: bool = True) -> None:
        """
        Recursively detect cycles in the subgraph starting from the current node.

        Parameters
        ----------
        current_node : Any
            The node currently being visited.
        visited : Dict[Any, bool]
            A dictionary keeping track of visited nodes.
        path : List[Any]
            A list to keep track of the current path.
        cycles : List[List[Any]]
            A list to accumulate all detected cycles.
        parent_node : Any, optional
            The parent node from which the current node was visited, by default None.
        is_undirected : bool, optional
            Indicates if the graph is undirected, by default True.

        Returns
        -------
        None
        """
        # Mark current node as visited
        visited[current_node] = True
        path.append(current_node)

        # Recursively visit all nodes along the branch from the current node
        for neighbour_node in self._nodes[current_node]['neighbours']:
            # Skip the edge leading back to the parent node (required for undirected graphs)
            if is_undirected and neighbour_node == parent_node:
                continue

            # If the neighbour/leaf node has not yet been visited, recursively visit its neighbour/leaf nodes
            if not visited[neighbour_node]:
                self._detect_cycles(neighbour_node, visited, path, cycles, current_node)

            # If current neighbour/leaf node has been visited in current branch path, there must be a cycle
            elif neighbour_node in path:
                cycle_start_index = path.index(neighbour_node)
                cycle = path[cycle_start_index:]
                if cycle not in cycles:
                    cycles.append(cycle)

        # After exploring all neighbours/leaf nodes in branch, current node is removed from path
        path.pop()

    def _dfs(self, current_node: Any, end_node: Any, predecessor: Dict[Any, Any], visited_nodes: set) -> bool:
        """
        Depth-first search (DFS) algorithm to recursively traverse graph from a given start node to given end node.

        Parameters
        ----------
        current_node : Any
            The current node for DFS.
        end_node : Any
            The target node for DFS.
        predecessor : Dict[Any, Any]
            A map of nodes to their predecessors.
        visited_nodes : set
            A set of visited nodes.

        Returns
        -------
        bool
            True if a path to end_node is found, else False.
        """
        visited_nodes.add(current_node)

        if current_node == end_node:
            return True

        for neighbour_node in self._nodes[current_node]['neighbours']:
            if neighbour_node not in visited_nodes:
                predecessor[neighbour_node] = current_node
                if self._dfs(neighbour_node, end_node, predecessor, visited_nodes):
                    return True

        return False

    def _dfs_connected_component(self, start_node: Any, visited_nodes: Set[Any]) -> None:
        """
        Perform a DFS traversal to find all nodes in the connected component nodes from a given start node.

        Parameters
        ----------
        start_node : Any
            The starting node for the DFS traversal.
        visited_nodes : Set[Any]
            The set of visited nodes to update.

        Returns
        -------
        None
        """
        stack = [start_node]
        while stack:
            current_node = stack.pop()
            if current_node not in visited_nodes:
                visited_nodes.add(current_node)
                for neighbour_node in self._nodes[current_node]['neighbours']:
                    if neighbour_node not in visited_nodes:
                        stack.append(neighbour_node)

    def _bfs(self, start_node: Any, end_node: Any) -> Union[List[Any], None]:
        """
        Breadth-first search (BFS) algorithm to traverse graph from a given starting node.

        Parameters
        ----------
        start_node : Any
            The starting node for BFS.
        end_node : Any
            The target node for BFS.

        Returns
        -------
        Union[List[Any], None]
            A list representing the shortest path if it exists, else None.
        """
        visited_nodes = set()
        predecessor = {start_node: None}

        self._fifo_queue.enqueue(start_node)
        visited_nodes.add(start_node)

        while not self._fifo_queue.empty():
            current_node = self._fifo_queue.get()

            if current_node == end_node:
                return self._reconstruct_path(predecessor, end_node)

            for neighbour_node in self._nodes[current_node]['neighbours']:
                if neighbour_node not in visited_nodes:
                    visited_nodes.add(neighbour_node)
                    self._fifo_queue.enqueue(neighbour_node)
                    predecessor[neighbour_node] = current_node

        return None

    @staticmethod
    def _reconstruct_path(predecessor: Dict[Any, Any], end_node: Any) -> List[Any]:
        """
        Reconstruct the shortest path from the predecessor map once end node has been found.

        Parameters
        ----------
        predecessor : Dict[Any, Any]
            A map of nodes to their predecessors.
        end_node : Any
            The ending node.

        Returns
        -------
        List
            The shortest path from start_node to end_node.
        """
        path = []
        current_node = end_node

        while current_node is not None:
            path.append(current_node)
            current_node = predecessor[current_node]

        # Reverse the path to get it from start to end
        path.reverse()
        return path

    @staticmethod
    def _has_traversable_neighbours(node: Any) -> bool:
        """
        Check if a node has traversable neighbors.

        Parameters
        ----------
        node : Any
            The node to check.

        Returns
        -------
        bool
            True if the node has neighbors, False otherwise.
        """
        return bool(node['neighbours'])


class UndirectedMolecularGraph(BaseMolecularGraph):
    """
    Undirected molecular graph data structure. Extends the BaseMolecularGraph parent class.

    Methods
    -------
    is_cyclic() -> bool
        Determine if the undirected molecular graph contains any cycles.
    """

    def is_cyclic(self) -> Union[List[List], None]:
        """
        Determine if the undirected molecular graph contains any cycles. If it does, return list of cycles, else return
        None

        Returns
        -------
        Union[List[List], None]
            Return List of cycles if any exist, else return None
        """
        # Initialise a `visited` dictionary to keep track of visited nodes. Initialise nodes as not visited
        visited: Dict = {node: False for node in self._nodes}
        path: List = []
        cycles: List = []

        # Iterate over all unvisited nodes and call _detect_cycles() helper method to recursively check for cycles
        for node in visited:
            # Only detect cycles if node is unvisited
            if not visited[node]:
                self._detect_cycles(node, visited, path, cycles)

        return cycles if cycles else None

    def _add_edge(self, from_node: Any, to_node: Any, **attr) -> None:
        """
        Add an edge between two nodes in the undirected graph. In undirected graph, edges are bidirectional, so the
        edge is added in both directions.

        Parameters
        ----------
        from_node : Any
            The starting node of the edge.
        to_node : Any
            The ending node of the edge.
        **attr : Dict
            Any additional edge attributes.

        Returns
        -------
        None

        Raises
        ------
        NonexistentNodeError
            If either the start node or the end node is not present in the graph.
        """
        # Validate that both `from_node` and `to_node` exist in the graph
        self._validate_nodes(from_node, to_node)

        self._nodes[from_node]['neighbours'].append(to_node)
        self._nodes[to_node]['neighbours'].append(from_node)
        self._edges[(from_node, to_node)] = attr


class DirectedMolecularGraph(BaseMolecularGraph):
    """
    Directed molecular graph data structure. Extends the BaseMolecularGraph parent class.

    Methods
    -------
    is_cyclic() -> bool
        Determine if the graph contains any cycles.
    topological_sort() -> List[Any]
        Perform topological sorting of the nodes in the graph.
    """

    def is_cyclic(self) -> Union[List[List], None]:
        """
        Determine if the directed molecular graph contains any cycles.

        Returns
        -------
        Union[List[List], None]
            Return List of cycles if any exist, else return None
        """
        # Initialise a `visited` dictionary to keep track of visited nodes. Initialise nodes as not visited
        visited: Dict = {node: False for node in self._nodes}
        path: List = []
        cycles: List = []

        # Iterate over all unvisited nodes and call _detect_cycles() helper method to recursively check for cycles
        for node in visited:
            # Only detect cycles if node is unvisited
            if not visited[node]:
                self._detect_cycles(node, visited, path, cycles, is_undirected=False)

        return cycles if cycles else None

    def topological_sort(self) -> List[Any]:
        """
        Perform topological sorting of the nodes in the graph.

        Returns
        -------
        list
            A list of nodes in topologically sorted order.

        Raises
        ------
        CycleError
            If the graph contains a cycle.
        """
        # A cyclic graph cannot have valid topological sorting. Raise error if graph contains cycles
        if self.is_cyclic():
            raise CycleError('A cyclic graph cannot have valid topological sorting')

        # Initialise a `visited` dictionary to keep track of visited nodes. Initialise nodes as not visited
        visited: Dict = {node: False for node in self._nodes}

        # Initialise stack that will be used to store the topologically sorted nodes
        stack: List = []

        # Iterate over all unvisited nodes and call _topological_sort_util() helper method to recursively store nodes
        # in topological sorted list
        for node in visited:
            # Only topologically sort unvisited nodes
            if not visited[node]:
                self._topological_sort_util(node, visited, stack)

        return stack

    def _add_edge(self, from_node: Any, to_node: Any, **attr) -> None:
        """
        Add an edge between two nodes in the directed graph. In directed graph, edges are unidirectional, so the
        edge is only added in the direction of the `to_node` node.

        Parameters
        ----------
        from_node : Any
            The starting node of the edge.
        to_node : Any
            The ending node of the edge.
        **attr : Dict
            Any additional node attributes.

        Returns
        -------
        None

        Raises
        ------
        NonexistentNodeError
            If either the start node or the end node is not present in the graph.
        """
        # Validate that both `from_node` and `to_node` exist in the graph
        self._validate_nodes(from_node, to_node)

        self._nodes[from_node]['neighbours'].append(to_node)
        self._edges[(from_node, to_node)] = attr

    def _topological_sort_util(self, current_node: Any, visited: Dict[Any, bool], stack: List) -> None:
        """
        Helper method to perform topological sorting.

        Parameters
        ----------
        current_node : Any
            The node currently being visited.
        visited : Dict[Any, bool]
            A dictionary keeping track of visited nodes.
        stack : List
            A list to store the topologically sorted nodes.


        Returns
        -------
        None
        """
        # Mark current node as visited
        visited[current_node] = True

        # Recursively visit all neighbour nodes along the branch from the current node
        for edge in self._nodes[current_node]:
            neighbour_node_tuple: Tuple = list(edge.items())[0]
            neighbour_node: Any = neighbour_node_tuple[0]

            if not visited[neighbour_node]:
                self._topological_sort_util(neighbour_node, visited, stack)

        # Once all neighbour nodes have been visited, insert node in current stack frame at the beginning of the
        # `stack` list
        stack.insert(0, current_node)
