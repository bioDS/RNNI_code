__author__ = 'Lena Collienne and @gavruskin'
# Methods for trees including writing nd reading from files and RNNI moves

import re

def read_trees_from_file():
    # Reads trees from tree_file and returns them as tree_list.
    # Trees assumed to be in tau-format.
    tree_file = input("What is the file with trees?\n")
    f = open(tree_file, 'r')
    tree_list = []
    for tree in f:
        tree_list.append(tree)
    f.close()
    return tree_list


def string_to_list_of_sets(s):
    # Converts a given string to a list of sets, if the string is representing a tree in our list representation
    s = s.rstrip() #Delete newline in the end of the string
    s = s.split('{')
    s.pop(0) # Because we split at '{' the first string does not contain information
    output = []
    for i in s:
        current_set = set()
        m = re.findall(r"(\d+)",i)
        for j in m:
            current_set.add(j)
        output.append(current_set)
    return output


def print_tree_list(tree_list):
    for tree in tree_list:
        print("list[%s] = %s" % (tree_list.index(tree), tree))


def get_num_taxa(tree):
        n = 0
        for i in tree:
            if len(i) > n:
                n = len(i)
        return n


def is_resolved(tree):
        if len(tree) == 2 * get_num_taxa(tree) - 1:
            return True
        elif len(tree) < 2 * get_num_taxa(tree) - 1:
            return False
        else:
            print("Not a tree--too many splits.")
            return


def is_tree(tree):
    # Returns True if tree is a tree and False otherwise.
    for i in range(len(tree)):
        for j in range(i+1,len(tree)):
            if tree[i] == tree[j]:
                print("Identical splits present.")
                return False
            elif tree[j].issubset(tree[i]):
                print("Refinement order is broken.")
                return False
            elif not tree[i].issubset(tree[j]) and len(tree[i].intersection(tree[j])) > 0:
                print("Splits with non-empty symmetric difference present.")
                return False
    return True


def sort_tree(tree):
    # Returns unique string representation of tree by sorting taxa and converting sets to trees.
    output = []
    for i in tree:
        j = sorted(list(i))
        output.append(set(j))
    return output


def swap_ranks(tree, i):
    # Swaps ranks of tree[i] and tree[i+1] if they are not adjacent by an edge.
    # Returns a list that consists of one tree.
    if i < len(tree) and len(tree[i].intersection(tree[i+1])) == 0:
        swapped_tree = list(tree)
        swapped_tree[i], swapped_tree[i+1] = swapped_tree[i+1], swapped_tree[i]
        return [swapped_tree]
    #else:
    #    print("Trying to swap unswappable.")


def make_nni(tree, i):
    # Returns a list consisting of two tree NNI adjacent to tree via edge (tree[i], tree[i+1]) if such edge exists.
    if i < len(tree) and tree[i].issubset(tree[i+1]):
        nni_tree1 = list(tree)
        nni_tree2 = list(tree)
        done = False
        j = i-1
        while not done:
            if tree[j].issubset(tree[i]):
                done = True
            else:
                j -= 1
        nni_tree1[i] = tree[j].union(tree[i+1].difference(tree[i]))
        nni_tree2[i] = tree[i].difference(tree[j]).union(tree[i+1].difference(tree[i]))
        return [nni_tree1, nni_tree2]
    #else:
    #   print("Trying to NNI inNNIable.")


def one_neighbourhood(tree):
    # Returns a list of trees that are one RNNI move apart from the given tree
    N = []
    nni = []
    for i in range(1,len(tree)-1):
        print(tree[i])
        rank = swap_ranks(tree, i)
        nni  = make_nni(tree,i)
        if rank != None:
            N = N + rank
        if nni!= None:
            N = N + nni
    return N

def one_neighbourhood_uRNNI(tree):
    # Returns a list of trees that are one RNNI move apart from the given tree, the given tree must be in RNNI but there are no moves performed on the first n sets (trivial cluster)
    N = []
    nni = []
    for i in range(int((len(tree)+1)/2), len(tree)-1):
        rank = swap_ranks(tree, i)
        nni  = make_nni(tree,i)
        if rank != None:
            N = N + rank
        if nni!= None:
            N = N + nni
    return N
