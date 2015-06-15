#!/usr/bin/env python3

import networkx as nx
import sys

def main():
    graph = nx.read_dot(sys.stdin)
    graph = graph.reverse()

    keep_nodes = {"all",
                  "raw/mcra-clones.all.fastq",
                  "raw/mcra.published.fn",
                  "unarchive.mk",
                  "mcra.published.fn",
                  "mcra-clones.all.fastq.mk"}
    drop_nodes = {".gitmodules", ".confirm-git-mangle", "venv"}
    out_degree = graph.out_degree()
    for node in graph.nodes_iter():
        path_parts = node.split("/")
        if node in keep_nodes:
            pass
        elif out_degree[node] == 0:
            drop_nodes.add(node)
        elif path_parts[0] in ("raw", "bin", ".git"):
            drop_nodes.add(node)
        else:
            pass
    graph.remove_nodes_from(drop_nodes)

    in_degree = graph.in_degree()
    for node in graph.node:
        if in_degree[node] == 0:
            graph.node[node]['shape'] = 'hexagon'
        rootdir = node.split("/")[0]
        if rootdir == 'raw':
            graph.node[node]['shape'] = 'box'
        elif rootdir == 'etc':
            pass
            # graph.node[node]['shape'] = 'circle'

    graph.graph = {}
    graph.graph['node'] = {'shape': 'plaintext'}
    # write_dot doesn't work, because it calls nx.to_agraph(graph).clear()
    # nx.write_dot(graph)
    a = nx.to_agraph(graph)
    a.write(sys.stdout)

if __name__ == "__main__":
    main()
