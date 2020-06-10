__author__ = 'Lena Collienne'
# Functions to compute the uRNNI graph (ultrametric Ranked Nearest Neighbour Interchange) as described in (Gavryushkin, Whidden, Matsen 2018) and computing shortest paths within this graph

import os.path
import os

from tree_generator import *
from graph import *


def uRNNI_graph_on_n_taxa(n):
    # Computes the uRNNI graph on n taxa
    uRNNI = Graph(set(), set())
    tree_list = [] #list containing all trees; in the end tree_list[i] is represented by vertex i in output graph
    all_trees_on_n_taxa(n)
    tree_file = open("all_trees_on_%s_taxa.txt" % n, 'r')
    # Example of tree on 3 taxa as given in file: [[1, 2], [1, 2, 3]]
    for line in tree_file:
        tree_list.append(string_to_tree(line.rstrip('\n')))
    tree_file.close()
    # Fill vertex set
    for i in range(0,len(tree_list)):
        uRNNI.vertices.add(i)
    for tree in tree_list:
        # Add all edges incident to tree
        nonu_tree = add_leaves(tree)
        # Add all edges incident to tree to uRNNI graph
        for neighbour in one_neighbourhood_uRNNI(nonu_tree):
            # First delete the sets containing just one leaf and then convert the tree to a string -> ultrametric tree u_neighbour
            u_neighbour = []
            for i in neighbour:
                if len(i)>1:
                    u_neighbour.append(i)
            # Add edge (tree, u_neighbour), if it is not already there
            edge1 = '(' + str(tree_list.index(tree)) + ',' + str(tree_list.index(u_neighbour)) + ')'
            edge2 = '(' +  str(tree_list.index(u_neighbour)) + ',' + str(tree_list.index(tree)) + ')'
            if edge1 not in uRNNI.edges and edge2 not in uRNNI.edges:
                uRNNI.edges.add(edge1)

    #print uRNNI graph in file
    try:
        os.makedirs("output/uRNNI_graphs")
    except OSError:
        if not os.path.isdir("output/uRNNI_graphs"):
            raise
    uRNNI.print_graph_to_file('output/uRNNI_graphs/uRNNI_edges_' + str(n) + '_taxa.txt')
    #print trees into file, number of line of tree equals vertex number that represents this tree in uRNNI graph
    f = open('output/uRNNI_graphs/tree_numbers_' + str(n) + '_taxa.txt','w')
    for tree in tree_list:
        f.write(tree_to_string(tree) + '\n')
    f.close()
    #output: uRNNI graph and tree_list
    return uRNNI, tree_list

def add_leaves(tree):
    # Add sets containing just one leaf to make ultrametric trees become non-ultrametric RNNI trees
    # needed to use RNNI moves as they are implemented for RNNI and not for uRNNI
    tree_with_leaves = []
    n = len(tree) + 1
    for j in range(1,n+1):
        tree_with_leaves.append(set([str(j)]))
    for k in tree:
        tree_with_leaves.append(set(k))
    return tree_with_leaves


def string_to_tree(S):
    # If a tree in our representation (list of lists) is converted to a string, this function converts it back to a tree
    string = S.replace('[' , '').split('],')
    tree = list(map(lambda s: s.split(', '),  string))
    for i in range(len(tree)):
        tree[i] = set(map(lambda x: int(x.strip(']]\n')),tree[i]))
    return tree


def tree_to_string(tree):
    # Convert a given tree as list of sets to list of lists saved as string
    # When saving lists of lists, the taxa in each list are sorted increasingly so that we compare two trees in this representation
    tree1 = list()
    for s in tree:
        l = list(s)
        l.sort()
        tree1.append(l)
    tree1 = str(tree1)
    return tree1

def shortest_uRNNI_path(start, end):
    # Computes a shortest path in uRNNI between the given trees start and end in our tree representation (list of sets) (by using shortest_path function for graphs)
    # Convert the given trees into the format we will use here
    start = convert_input_to_tree(start)
    end = convert_input_to_tree(end)
    n = len(start)+1
    # Compute uRNNI graph on n taxa, if it doesn't exist already
    if os.path.exists('output/uRNNI_graphs/tree_numbers_%s_taxa.txt' %n):
        # print("Start reading uRNNI graph on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
        uRNNI = read_uRNNI_graph('output/uRNNI_graphs/tree_numbers_%s_taxa.txt' %n, 'output/uRNNI_graphs/uRNNI_edges_%s_taxa.txt' %n)
        # print("Done reading uRNNI graph on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
    else:
        # Compute the uRNNI graph on n taxa
        print("Start generating uRNNI graph on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
        uRNNI = uRNNI_graph_on_n_taxa(n)
        print("Done generating uRNNI graph on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
    # Get integers representing the given start and end tree in the RNNI graph
    start_node = uRNNI[1].index(str(start))
    end_node = uRNNI[1].index(str(end))
    # Compute and return distance and a shortest path between start and end tree using the function 'shortest_path' of the class 'graph'
    sp = uRNNI[0].shortest_path(start_node,end_node)
    distance = sp[0]
    path = sp[1]
    trees_on_path = []
    for i in path:
        trees_on_path.append(uRNNI[1][i])
    return (distance, trees_on_path)

def convert_input_to_tree(tree):
    # Converts the given tree representation (list of sets) into the one we will use in the program (list of lists with increasingly ordered numbers)
    tree_list = []
    current_list = []
    for i in tree:
        current_list = []
        for j in i:
            current_list.append(int(j))
        tree_list.append(sorted(current_list))
    return tree_list

def read_uRNNI_graph(nodefile, edgefile):
    #nodefile: trees saved as list of lists, edgefile: edges saved as (v,w), v and w are the integers corresponding to the trees in edgefile
    tree_list = []
    nfile = open(nodefile, 'r')
    i = 0
    vertices = set() #vertex set of the uRNNI graph, to be filled with integers 0,...,N:= number of trees on n taxa
    # First read vertices
    for line in nfile:
        vertices.add(i)
        tree_list.append(string_to_tree(line.rstrip('\n')))
        i += 1
    nfile.close()
    # Second read edges
    edges = set()
    efile = open(edgefile, 'r')
    for line in efile:
        edges.add(line.rstrip('\n'))
    uRNNI = [Graph(vertices, edges), tree_list]
    return uRNNI
