from __future__ import annotations

from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd


def load_string_graph(links_path: str | Path, minimum_score: int = 700) -> nx.Graph:
    frame = pd.read_csv(links_path, sep=r"\s+")
    required = {"protein1", "protein2", "combined_score"}
    if missing := required.difference(frame.columns):
        raise ValueError(f"STRING links file is missing {sorted(missing)}")
    frame = frame.loc[frame["combined_score"] >= minimum_score]
    graph = nx.from_pandas_edgelist(frame, "protein1", "protein2", edge_attr="combined_score")
    graph.remove_edges_from(nx.selfloop_edges(graph))
    return graph


def induced_subgraph_statistics(graph: nx.Graph, proteins: set[str]) -> dict[str, float]:
    subgraph = graph.subgraph(proteins).copy()
    return {
        "n_input_proteins": len(proteins),
        "n_mapped_proteins": sum(node in graph for node in proteins),
        "n_nodes": subgraph.number_of_nodes(),
        "n_edges": subgraph.number_of_edges(),
        "density": nx.density(subgraph) if subgraph.number_of_nodes() > 1 else 0.0,
        "largest_component": max((len(component) for component in nx.connected_components(subgraph)), default=0),
    }


def degree_matched_edge_test(graph: nx.Graph, proteins: set[str], permutations: int = 1000, seed: int = 20260814) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    targets = [node for node in proteins if node in graph]
    observed = graph.subgraph(targets).number_of_edges()
    degree = dict(graph.degree())
    nodes = np.array(list(graph.nodes()))
    node_degrees = np.array([degree[node] for node in nodes])
    bins = np.quantile(node_degrees, np.linspace(0, 1, 11))
    target_bins = np.digitize([degree[node] for node in targets], bins[1:-1], right=True)
    null = []
    for _ in range(permutations):
        sampled = []
        for bin_id in target_bins:
            pool = nodes[np.digitize(node_degrees, bins[1:-1], right=True) == bin_id]
            sampled.append(rng.choice(pool))
        null.append(graph.subgraph(set(sampled)).number_of_edges())
    null = np.asarray(null)
    return {
        "observed_edges": float(observed),
        "null_mean_edges": float(null.mean()),
        "empirical_p_value": float((1 + np.sum(null >= observed)) / (permutations + 1)),
        "permutations": permutations,
    }

