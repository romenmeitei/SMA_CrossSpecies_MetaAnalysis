import networkx as nx

from sma_xspecies.network import induced_subgraph_statistics


def test_induced_network_keeps_edges_between_arbitrary_input_partitions():
    graph = nx.Graph()
    graph.add_edges_from([("A", "B"), ("B", "C"), ("C", "D")])
    stats = induced_subgraph_statistics(graph, {"A", "B", "C", "D"})
    assert stats["n_edges"] == 3
    assert stats["largest_component"] == 4

