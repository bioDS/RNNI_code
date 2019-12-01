__author__ = "Lena Collienne"

#functions on caterpillar trees, specifically computing distances and shortest paths between caterpillar trees

from tree_generator import *
import copy

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


def caterpillar_list(tree):
    # Returns a list representing a caterpillar tree as introduced in Caterpillar Sort (Section  1.2)
    n = get_num_taxa(tree)
    j = 1
    clist = []
    clist.append(int(min(tree[0])))
    clist.append(int(max(tree[0])))
    # Add all remaining taxa to the list
    for i in range(1,len(tree)):
        clist.append(int((tree[i].difference(tree[i-1])).pop()))
    return clist


def caterpillar_list_to_tree(clist):
    tree = []
    for i in range(len(clist)-1): #Add cherry taxa to all sets
        tree.append(set([clist[0], clist[1]]))
    for i in range(2,len(clist)):
        for j in range(i-1, len(tree)):
            tree[j].add(clist[i])
    return(tree)


def shortest_caterpillar_path(tree1, tree2):
    # Returns a shortest caterpillar path between two caterpillar trees (ultrametric)
    if is_caterpillar(tree1) and is_caterpillar(tree2) and get_num_taxa(tree1) == get_num_taxa(tree2):
        clist1 = caterpillar_list(tree1)
        clist2 = caterpillar_list(tree2)
        current_tree = clist1
        n = len(clist1)
        cpath = [copy.copy(current_tree)] # computed path with trees in caterpillar list representation
        for j in range(n-1, 1, -1):
            for i in range (0, j):
                if current_tree[i] == clist2[j]:
                    x = current_tree[i]
                    current_tree[i] = current_tree[i+1]
                    current_tree[i+1] = x
                    cpath.append(copy.copy(current_tree))
        path = [] # Output path with trees in cluster representation
        for i in cpath:
            path.append(caterpillar_list_to_tree(i))
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
