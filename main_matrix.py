from tree import *
from graph import *

__author__ = '@gavruskin'


def delta(tree):
    # Returns 1-neighbourhood of tree as a list, delta(tree)[0] == tree.
    # tree is assumed to be a list of splits, which are sets of taxa.
    output = [tree]
    # The following adds trees obtained from tree by swapping ranks of nodes:
    for i in range(len(tree) - 1):
        if len(tree[i].intersection(tree[i+1])) == 0:
            output += swap_ranks(tree, i)
        elif tree[i].issubset(tree[i+1]):
            output += make_nni(tree, i)
    return output


def graph_on_trees(tree_list):
    # Returns tau-graph on trees from tree_list.
    # tree_list must contain no repetitions
    graph = Graph(set(), [])
    for i in range(len(tree_list)):
        graph.vertices.add(i)
        for tree in delta(tree_list[i]):
            for j in range(i):
                if sort_tree(tree_list[j]) == sort_tree(tree):
                    graph.edges.append(set([i,j]))
    return graph
