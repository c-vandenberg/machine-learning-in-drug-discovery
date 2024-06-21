from typing import Any, Union, List, Set, Dict, Tuple
from queues import FifoQueue
from helpers.exception import CycleError


class BaseMolecularGraph:
    """
   Base class for a molecular graph data structure.

   Methods
   -------
   nodes() -> Dict[Union[int, float], List[Dict]]
       Getter for the graph nodes and their edges.
   add_node(node: Any)
       Add node to graph if it doesn't already exist.
   add_edge(from_node: Any, to_node: Any, weight: Union[int, float, None] = None)
       Add edge between two nodes in graph.
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
        self._visited_nodes = set()
        self._traversal_order = []
        self._has_path = False

    @property
    def nodes(self) -> dict:
        """
        Getter for the graph nodes

        Returns
        -------
        dict
            A dictionary where keys are nodes and the values are sub dictionaries of the node attributes
        """
        return self._nodes

    @property
    def edges(self) -> dict:
        """
        Getter for the graph edges

        Returns
        -------
        dict
            A dictionary where keys are tuples representing the graph edges, and the values are sub dictionaries of the
            edge attributes
        """
        return self._edges

    def add_node(self, node: Any, **attr):
        """
        Add node to graph if it doesn't already exist.

        Parameters
        ----------
        node : Any
            The node to be added to the graph.
        **attr : dict
            Any additional node attributes.
        """
        if node not in self._nodes:
            # Nodes are defined as a dictionary with 'neighbors' key and attributes
            self._nodes[node] = {'neighbours': [], **attr}
        else:
            # If node already exists, update the attributes
            self._nodes[node].update(attr)

    def add_edge(self, from_node: Any, to_node: Any, **attr):
        """
        Add edge between two nodes in graph.

        Parameters
        ----------
        from_node : Any
            The starting node of the edge.
        to_node : Any
            The ending node of the edge.
        **attr : dict
            Any additional edge attributes.
        """
        self._add_edge(from_node, to_node, **attr)

    def _add_edge(self, from_node: Any, to_node: Any, **attr):
        """
        Add edge between two nodes in graph (to be implemented/overridden by subclasses).

        Parameters
        ----------
        from_node : Any
            The starting node of the edge.
        to_node : Any
            The ending node of the edge.
        **attr : dict
            Any additional edge attributes.
        """
        pass

    def _validate_nodes(self, start_node: Any, end_node: Union[None, Any] = None):
        """
        Validate that the nodes exist in the graph.

        Parameters
        ----------
        start_node : Any
            The starting node.
        end_node : Union[None, Any]. Optional
            The ending node, by default None.

        Raises
        ------
        ValueError
            If either the start node or the end node is not present in the graph.
        """
        if start_node not in self._nodes or (end_node not in self._nodes and end_node):
            missing_node: Any = start_node if start_node not in self._nodes else end_node

            raise ValueError(f'Node {missing_node} not present in Graph')

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
        list or None
            A list representing the path if it exists, otherwise None.
        """
        self._validate_nodes(start_node, end_node)
        path: Union[List, None] = None
        self._dfs(start_node, end_node)

        if self._has_path:
            path = self._traversal_order[:]
            self._has_path = False

        self._visited_nodes.clear()
        self._traversal_order.clear()

        return path

    def _dfs(self, node: Any, end_node: Union[None, Any] = None) -> None:
        """
        Depth-first search (DFS) algorithm to traverse graph from a given starting node.

        Parameters
        ----------
        node : Any
            The starting node for DFS.
        end_node : Union[None, Any]. Optional
            The target node for DFS, by default None.
        """
        self._visited_nodes.add(node)
        self._traversal_order.append(node)

        if node == end_node:
            self._has_path = True
            return

        for neighbour_node in self._nodes[node]['neighbours']:
            if neighbour_node not in self._visited_nodes:
                self._dfs(neighbour_node, end_node)

    def connected_components(self) -> Set:
        """
        Return all connected components in the graph using depth-first search (DFS) algorithm to traverse graph.

        Returns
        -------
        set
            A set of nodes that have neighbors (connected components).
        """
        for node in self._nodes:
            if self._has_traversable_neighbours(self._nodes[node]) and node not in self._visited_nodes:
                self._dfs(node)

        connected_components: Set = set(self._visited_nodes)
        self._visited_nodes.clear()
        self._traversal_order.clear()

        return connected_components

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
        list or None
            A list representing the shortest path if it exists, otherwise None.
        """
        self._validate_nodes(start_node, end_node)
        path: Union[List, None] = None
        self._bfs(start_node, end_node)

        if self._has_path:
            path = self._traversal_order[:]
            self._has_path = False

        self._visited_nodes.clear()
        self._traversal_order.clear()

        return path

    def _bfs(self, start_node: Any, end_node: Union[None, Any] = None) -> None:
        """
        Breadth-first search (BFS) algorithm to traverse graph from a given starting node.

        Parameters
        ----------
        start_node : Any
            The starting node for BFS.
        end_node : Union[None, Any]. Optional
            The target node for BFS, by default None.
        """
        if len(self._visited_nodes) != 0:
            self._visited_nodes.clear()

        if len(self._traversal_order) != 0:
            self._traversal_order.clear()

        self._fifo_queue.enqueue(start_node)
        self._visited_nodes.add(start_node)
        self._traversal_order.append(start_node)

        while self._fifo_queue:
            current_node = self._fifo_queue.dequeue()

            if current_node == end_node:
                self._has_path = True
                return

            for edge in self._nodes[current_node]:
                neighbour_node: Any = list(edge.keys())[0]
                if neighbour_node not in self._visited_nodes:
                    self._fifo_queue.enqueue(neighbour_node)
                    self._visited_nodes.add(neighbour_node)
                    self._traversal_order.append(neighbour_node)

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
        return True if len(node) != 0 else False

    def __str__(self):
        """
        Return the string representation of the graph.

        Returns
        -------
        str
            A string representation of the nodes and their edges.
        """
        return str(self._nodes)


class UndirectedMolecularGraph(BaseMolecularGraph):
    """
    Undirected molecular graph data structure. Extends the BaseMolecularGraph parent class.

    Methods
    -------
    is_cyclic() -> bool
        Determine if the graph contains any cycles.
    """

    def _add_edge(self, from_node: Any, to_node: Any, weight: Union[int, float, None] = None, **attr):
        """
        Add an edge between two nodes in the undirected graph. In undirected graph, edges are bidirectional, so the
        edge is added in both directions.

        Parameters
        ----------
        from_node : Any
            The starting node of the edge.
        to_node : Any
            The ending node of the edge.
        **attr : dict
            Any additional edge attributes.
        """
        # Validate that both `from_node` and `to_node` exist in the graph
        self._validate_nodes(from_node, to_node)

        self._nodes[from_node]['neighbours'].append(to_node)
        self._nodes[to_node]['neighbours'].append(from_node)
        self._edges[(from_node, to_node)] = attr

    def is_cyclic(self) -> Union[List, None]:
        """
        Determine if the graph contains any cycles. If it does, return list of cycles, else return None

        Returns
        -------
        Union[List, None]
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

    def _detect_cycles(self, current_node: Any, visited: Dict[Any, bool], path: List, cycles: List,
                       parent_node: Any = None):
        """
        Recursively detect cycles in the subgraph starting from the current node.

        Parameters
        ----------
        current_node : Any
            The node currently being visited.
        visited : Dict[Any, bool]
            A dictionary keeping track of visited nodes.
        parent_node : Any, optional
            The parent node from which the current node was visited, by default None.

        Returns
        -------
        bool
            True if a cycle is detected, False otherwise.
        """
        # Mark current node as visited
        visited[current_node] = True
        path.append(current_node)

        # Recursively visit all nodes along the branch from the current node
        for neighbour_node in self._nodes[current_node]['neighbours']:
            # Skip the edge leading back to the parent node (required for undirected graphs)
            if neighbour_node == parent_node:
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

    def _add_edge(self, from_node: Any, to_node: Any, weight: Union[int, float, None] =  None):
        """
        Add an edge between two nodes in the directed graph. In undirected graph, edges are unidirectional, so the
        edge is only added in the direction of the `to_node` node.

        Parameters
        ----------
        from_node : Any
            The starting node of the edge.
        to_node : Any
            The ending node of the edge.
        weight : Union[int, float, None]
            The weight of the edge. Optional
        """
        # Validate that both `from_node` and `to_node` exist in the graph
        self._validate_nodes(from_node, to_node)

        self._nodes[from_node].append(
            {to_node: weight}
        )

    def is_cyclic(self):
        """
        Determine if the graph contains any cycles.

        Returns
        -------
        bool
            True if the graph contains a cycle, False otherwise.
        """
        # Initialise a `visited` dictionary to keep track of visited nodes. Initialise nodes as not visited
        visited: Dict = {node: False for node in self._nodes}

        # Initialise a `recursion_stack` dictionary to keep track of visited nodes in current recursion stack.
        # Initialise nodes as not visited
        recursion_stack = {node: False for node in self._nodes}

        # Iterate over all unvisited nodes and call _detect_cycles() helper method to recursively check for cycles
        for node in visited:
            # Only detect cycles if node is unvisited
            if not visited[node]:
                if self._detect_cycles(node, visited, recursion_stack):
                    return True

        return False

    def _detect_cycles(self, current_node: Any, visited: Dict[Any, bool], recursion_stack: Union[Dict, bool]):
        """
        Recursively detect cycles in the subgraph starting from the current node.

        Parameters
        ----------
        current_node : Any
            The node currently being visited.
        visited : Dict[Any, bool]
            A dictionary keeping track of visited nodes.
        recursion_stack : Union[Dict, bool]
            A dictionary keeping track of the recursion stack.

        Returns
        -------
        bool
            True if a cycle is detected, False otherwise.
        """
        # Mark current node as visited both overall and within to recursion stack
        visited[current_node] = True
        recursion_stack[current_node] = True

        # Recursively visit all nodes along the branch from the current node
        for edge in self._nodes[current_node]:
            neighbour_node_tuple: Tuple = list(edge.items())[0]
            neighbour_node: Any = neighbour_node_tuple[0]

            # If the neighbouring node has not yet been visited, recursively visit its neighbouring nodes
            if not visited[neighbour_node]:
                if self._detect_cycles(neighbour_node, visited, recursion_stack):
                    return True

            # If current neighbour is in recursive stack then we have already visited it in this branch search and so
            # there must be a cycle in the graph
            elif recursion_stack[neighbour_node]:
                return True

        recursion_stack[current_node] = False
        return False

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

    def _topological_sort_util(self, current_node: Any, visited: Dict[Any, bool], stack: List):
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
