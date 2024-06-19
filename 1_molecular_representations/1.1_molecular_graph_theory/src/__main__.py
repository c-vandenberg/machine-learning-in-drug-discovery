#!/usr/bin/env python3

from queues import FifoQueue
from graphs import DirectedGraph, UndirectedGraph
from typing import Union, List, Set


def main():
    # Create directed graph
    directed_fifo_queue: FifoQueue = FifoQueue()
    my_directed_graph: DirectedGraph = DirectedGraph(directed_fifo_queue)

    # Add nodes
    my_directed_graph.add_node("A")
    my_directed_graph.add_node("B")
    my_directed_graph.add_node("C")
    my_directed_graph.add_node("D")
    my_directed_graph.add_node("E")
    my_directed_graph.add_node("F")
    my_directed_graph.add_node("G")
    my_directed_graph.add_node("H")
    my_directed_graph.add_node("I")

    # Add edges
    my_directed_graph.add_edge("A", "B")
    my_directed_graph.add_edge("B", "C")
    my_directed_graph.add_edge("B", "G")
    my_directed_graph.add_edge("C", "G")
    my_directed_graph.add_edge("G", "A")

    directed_graph_connected_nodes: Set = my_directed_graph.connected_components()
    directed_has_cycles: bool = my_directed_graph.is_cyclic()

    directed_graph_shortest_a_g_path: Union[List, None] = my_directed_graph.find_shortest_path(
        "A", "G"
    )

    # Create undirected graph
    undirected_fifo_queue: FifoQueue = FifoQueue()
    my_undirected_graph = UndirectedGraph(undirected_fifo_queue)

    # Add nodes
    my_undirected_graph.add_node("A")
    my_undirected_graph.add_node("B")
    my_undirected_graph.add_node("C")
    my_undirected_graph.add_node("D")
    my_undirected_graph.add_node("E")
    my_undirected_graph.add_node("F")
    my_undirected_graph.add_node("G")
    my_undirected_graph.add_node("H")
    my_undirected_graph.add_node("I")

    # Add edges
    my_undirected_graph.add_edge("A", "B")
    my_undirected_graph.add_edge("B", "C")
    my_undirected_graph.add_edge("C", "D")
    my_undirected_graph.add_edge("G", "A")

    undirected_graph_a_c_path: Union[List, None] = my_undirected_graph.find_path("A", "C")
    undirected_graph_a_g_path: Union[List, None] = my_undirected_graph.find_path("A", "G")
    undirected_has_cycles: bool = my_undirected_graph.is_cyclic()


if __name__ == '__main__':
    main()
