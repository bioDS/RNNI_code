import os.path
import random
import time
from tree import *


__author__ = '@gavruskin'


def is_tree_index(tree_index):
    # tree_index is a sequence of (ordered) pairs of indices that defines
    # the order in which nodes of the tree coalesce.
    # Important: tree_index[i][0] < tree_index[i][1]
    # Returns True if tree_index is a tuple that can be converted into tree,
    # False otherwise. Note that len(tree_index) must be taxa_number - 1.
    # NOT TESTED
    n = len(tree_index) + 1
    for i in range(len(tree_index)):
        if tree_index[i][0] > n - i or tree_index[i][1] > n - i:
            return False
    return True


def index_to_tree(tree_index):
    # Returns the tree on taxa 1,...,n, where n = len(tree_index) + 1, of index tree_index.
    running_index = []
    output_tree = []
    for i in range(len(tree_index) + 1):
        running_index.append([i+1])
    for i in range(len(tree_index)):
        left_index = tree_index[i][0]
        right_index = tree_index[i][1]
        left = running_index[left_index]
        right = running_index[right_index]
        running_index.append(left + right)
        output_tree.append(set(sorted(left + right)))
        del running_index[left_index]
        del running_index[right_index - 1]
    return output_tree


def is_final(index, n):
    # True if index of the tree on n taxa is final.
    for i in range(n-1):
        pair = index[i]
        if pair[1] != n - i or pair[0] != n - i - 1:
            return False
    return True


def next_tree_index(index, n):
    # Modifies the index of the tree on n taxa following the tree with the given index.
    # If the index is final, returns 0.
    if is_final(index, n):
        return 0
    for i in range(n-1):
        pair = index[i]
        if pair[1] < n - i:
            index[i] = [index[i][0], index[i][1] + 1]
            for j in range(i):
                index[j] = [1, 2]
            return index
        elif pair[0] < n - i - 1:
            index[i] = [index[i][0] + 1, index[i][0] + 2]
            for j in range(i):
                index[j] = [1, 2]
            return index
        elif pair[1] - pair[0] > 1:
            index[i] = [index[i][0] + 1, index[i][1]]
            for j in range(i):
                index[j] = [1, 2]
            return index


def all_trees_on_n_taxa(n):
    # Writes all trees on n taxa one by one into all_trees_1by1_as_indices.txt
    if not os.path.isfile("./all_trees_as_indices_on_%s_taxa.txt" % n):
        index_file = open("all_trees_as_indices_on_%s_taxa.txt" % n, "w")
        first_tree_index = []
        for i in range(n-1):
            first_tree_index.append([1, 2])
        index_file.write(str(first_tree_index) + "\n")
    else:
        print("The index-file is not empty")
        return
    if os.path.isfile("./all_trees_on_%s_taxa.txt" % n):
        print("The tree-file is not empty")
        return
    print("Start generating trees on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
    tree_file = open("all_trees_on_%s_taxa.txt" % n, "w")
    tree = index_to_tree(first_tree_index)
    tree_file.write(str(tree) + "\n")
    last_tree_index = first_tree_index
    while not is_final(last_tree_index, n):
        new_tree_index = next_tree_index(last_tree_index, n)
        index_file.write(str(new_tree_index) + "\n")
        tree = index_to_tree(new_tree_index)
        tree_file.write(str(tree) + "\n")
        last_tree_index = new_tree_index
    index_file.close()
    tree_file.close()
    print("Trees on %s taxa are ready on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))


def generate_random_trees(n, k):
    # Generates a random list of k trees on n taxa, repetitions possible.
    output = []
    for i in range(k):
        new_index = []
        for j in range(n-1):
            a = random.randint(0, n - 1 - j - 1)
            b = random.randint(a + 1, n - 1 - j)
            new_index.append([a, b])
        output.append(index_to_tree(new_index))
    return output


def drop_repetitions(tree_list):
    for tree in tree_list:
        tree = sort_tree(tree)
    output = []
    for i in range(len(tree_list)):
        is_seen = False
        for j in range(i):
            if tree_list[j] == tree_list[i]:
                is_seen = True
        if not is_seen:
            output.append(tree)
    return output
