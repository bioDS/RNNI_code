__author__ = "Lena Collienne"

#functions on caterpillar trees, specifically computing distances and shortest paths between caterpillar trees

from tree_generator import *


def is_caterpillar(tree):
    # Returns True if the given tree is a caterpillar tree and False otherwise
    if is_tree(tree) == False:
        return False
    else:
        k = -1
        for i in range(len(tree)-1):
            if len(tree[i]) > 1:
                k = i
            if len(tree[i+1]) > 1 and k != -1 and tree[k].intersection(tree[i+1]) != tree[k]:
                return False
    return True



def rank_list(tree):
    # Returns a list containing for each taxon its rank and the rank of its parent (at position 'P' + taxon)
    n = get_num_taxa(tree)
    j = 1
    rank = {}
    for i in range(1, n+1):
        rank[str(i)] = 0
        rank['P' + str(i)] = 0
    for i in range (len(tree)):
        if len(tree[i]) == 1:
            for x in tree[i]:
                rank[x] = j
        if len(tree[i]) > 1:
            for x in tree[i]:
                if rank['P' + x] == 0:
                    rank['P' + x] = j
        j = j + 1
    return rank



def shortest_caterpillar_path(tree1, tree2):
    # Returns a shortest caterpillar path between two caterpillar trees if all taxa have ranks smaller than all internal nodes
    path = [tree1]
    rank2 = rank_list(tree2)
    if is_caterpillar(tree1) and is_caterpillar(tree2) and get_num_taxa(tree1) == get_num_taxa(tree2):
        current_tree = tree1
        for j in range(len(current_tree), 1, -1):
            for i in range (0, j-1):
                # Special case: taxa
                if len(current_tree[i]) == 1:
                    taxon1 = current_tree[i].difference({}).pop()
                    # Swap ranks of two taxa if necessary
                    if len(current_tree[i+1]) == 1:
                        taxon2 = current_tree[i+1].difference({}).pop()
                        if (rank2[taxon1] > rank2[taxon2]) and swap_ranks(current_tree, i) != None:
                            current_tree = swap_ranks(current_tree, i)[0]
                            path = path + [current_tree]

                    # Swap ranks of taxon and internal node with higher rank (not parent)
                    elif taxon1 not in current_tree[i+1]:
                        # Find taxon which is child of current_tree[i+1]
                        k = i-1
                        while len(current_tree[k].intersection(current_tree[i+1])) == 0:
                            k = k - 1
                        taxon2 = current_tree[i+1].difference(current_tree[k]).pop()
                        if rank2[taxon1] > rank2['P' + taxon2]:
                            current_tree = swap_ranks(current_tree, i)[0]
                            path = path + [current_tree]

                # Rank swap of internal node and taxon where taxon has the higher rank
                elif len(current_tree[i]) > 1 and len(current_tree[i+1]) == 1:
                    taxon1 = current_tree[i+1].difference({}).pop()
                    k = i + 1
                    while (current_tree[k] == 1):
                        k = k + 1
                    taxon2 = current_tree[i].difference(current_tree[k]).pop()
                    if (rank2['P' + taxon2] > rank2[taxon1]):
                        current_tree = swap_ranks(current_tree, i)[0]
                        path = path + [current_tree]

                # Special case: cherry - two taxa with equal ranks
                elif len(current_tree[i]) == 2:
                    cherry = current_tree[i].difference({})
                    taxon1 = cherry.pop()
                    taxon2 = cherry.pop()
                    taxon3 = current_tree[i+1].difference(current_tree[i]).pop()
                    aux_trees = make_nni(current_tree, i)
                    aux_rank = rank_list(aux_trees[0])
                    if rank2['P' + taxon1] > rank2['P' + taxon3] and rank2['P' + taxon1] > rank2['P' + taxon2] and make_nni(current_tree, i) != None:
                        if (aux_rank['P' + taxon1] > aux_rank['P' + taxon2]):
                            current_tree = aux_trees[0]
                        else:
                            current_tree = aux_trees[1]
                        path = path + [current_tree]
                    elif rank2['P' + taxon2] > rank2['P' + taxon3] and rank2['P' + taxon1] < rank2['P' + taxon2] and make_nni(current_tree, i) != None:
                        if (aux_rank['P' + taxon1] < aux_rank['P' + taxon2]):
                            current_tree = aux_trees[0]
                        else:
                            current_tree = aux_trees[1]
                        path = path + [current_tree]

                # Regular case: nni move on internal nodes
                elif len(current_tree[i]) > 2:
                    taxon1 = current_tree[i].difference(current_tree[i-1]).pop()
                    taxon2 = current_tree[i+1].difference(current_tree[i]).pop()
                    if rank2['P' + taxon1] > rank2['P' + taxon2] and make_nni(current_tree, i) != None:
                        current_tree = make_nni(current_tree, i)[0]
                        path = path + [current_tree]
        return path
    else:
        print('Input trees are not caterpillars with equal number of taxa')



def caterpillar_distance(tree1, tree2):
    # Returns the distance between two caterpillar trees by using shortest_caterpillar_path
    if is_caterpillar(tree1) and is_caterpillar(tree2) and get_num_taxa(tree1) == get_num_taxa(tree2):
        path = shortest_caterpillar_path(tree1, tree2)
        return (len(path)-1)
    else:
        print('Input trees are not caterpillars with equal number of taxa')
