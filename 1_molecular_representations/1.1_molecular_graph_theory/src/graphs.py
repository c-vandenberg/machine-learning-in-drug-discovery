from queue import Queue, LifoQueue, PriorityQueue
from collections import defaultdict
from typing import Any, Union, List, Set, Dict, Tuple
from queues import FifoQueue
from helpers.exception import CycleError


class BaseGraph:
    """
   A base class for representing a graph structure.

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
        self._nodes = defaultdict(list)
        self._visited_nodes = set()
        self._traversal_order = []
        self._has_path = False

    @property
    def nodes(self) -> Dict[Union[int, float], List[Dict]]:
        """
        Getter for the graph nodes and their edges.

        Returns
        -------
        dict
            A dictionary where keys are nodes and values are lists of dictionaries representing edges.
        """
        return self._nodes

    def add_node(self, node: Any):
        """
        Add node to graph if it doesn't already exist.

        Parameters
        ----------
        node : Any
            The node to be added to the graph.
        """
        if node not in self._nodes:
            # Nodes are defined as dictionary key with the edges being a list of neighbouring nodes
            self._nodes[node] = []

    def add_edge(self, from_node: Any, to_node: Any, weight: Union[int, float, None] = None):
        """
        Add edge between two nodes in graph.

        Parameters
        ----------
        from_node : Any
            The starting node of the edge.
        to_node : Any
            The ending node of the edge.
        weight : Union[int, float, None]. Optional
            The weight of the edge, by default None.
        """
        self._add_edge(from_node, to_node, weight)

    def _add_edge(self, from_node: Any, to_node: Any, weight: Union[int, float, None] = None):
        """
        Add edge between two nodes in graph (to be implemented/overridden by subclasses).

        Parameters
        ----------
        from_node : Any
            The starting node of the edge.
        to_node : Any
            The ending node of the edge.
        weight : Union[int, float, None]
            The weight of the edge. Optional
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

        for edge in self._nodes[node]:
            neighbour_node: Any = list(edge.keys())[0]
            if neighbour_node not in self._visited_nodes:
                self._dfs(neighbour_node, end_node)

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

    """ Use BFS algorithm to find the shortest find path between two nodes """
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

    def connected_components(self) -> Set:
        """
        Return all connected components in the graph.

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
