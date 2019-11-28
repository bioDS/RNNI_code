__author__ = 'Lena Collienne'
# Implementation of algorithm FINDPATH in paper

from uRNNI_graph import *

def findpath(tree1, tree2):
    # Returns a path between tree1 and tree2 (cluster representation)
    # Can be used as approximation for the RNNI distance
    n = len(tree1) + 1
    path = [] # output path between tree1 and tree2
    tree1 = add_leaves(tree1)
    tree2 = add_leaves(tree2)
    distance = 0
    current_tree = tree1
    path.append(current_tree)
    for i in range(n,len(tree1)-1):
        # Find mrca(tree2[i]) in current tree (with rank j)
        j = i
        while current_tree[j].issuperset(tree2[i]) == False:
            j += 1
        # Move taxa down to position in end tree
        for k in range (j-1, i-1, -1):
            # With a NNI move
            if current_tree[k].issubset(current_tree[k+1]):
                nni_list = make_nni(current_tree, k)
                if nni_list[0][k].issuperset(tree2[i]):
                    current_tree = nni_list[0]
                else:
                    current_tree = nni_list[1]
            # With a rank move
            else:
                current_tree = swap_ranks(current_tree, k)[0]
            distance = distance + 1
            path.append(current_tree)
    return path
