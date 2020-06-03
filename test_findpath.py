__author__ = 'Lena Collienne, Kieran Elmes'

#checks whether FINDPATH computes correct distances
#works on small trees with up to 6 taxa

from findpath import *
from os.path import exists
from itertools import repeat
from multiprocessing import Pool
from multiprocessing import sharedctypes
from tqdm import tqdm
from p_tqdm import p_map
import re
import pickle
import os

D = np.zeros((1,1))
l = 0
uRNNI = 0
def check_tree(tree1):
    global D
    global l
    global uRNNI
    for tree2 in uRNNI[1]:
        findpath_len = len(findpath(string_to_tree(str(tree1)),string_to_tree(str(tree2)))) - 1
        actual_len = D[l[uRNNI[1].index(tree1)]][l[uRNNI[1].index(tree2)]]
        if (findpath_len != actual_len):
            print("{} != {}, Failed on:".format(findpath_len, actual_len))
            print(tree1, tree2)
            return False
    return True

def check_findpath(n):
    global D
    global l
    global uRNNI
    try:
        os.makedirs("output/distance_matrices")
    except OSError:
        if not os.path.isdir("output/distance_matrices"):
            raise
    # Check if FINDPATH computes exact distances by comparing these distances with exact ones (Seidel on RNNI graph)
    if exists('output/uRNNI_graphs/tree_numbers_%s_taxa.txt' %n):
        uRNNI = read_uRNNI_graph('output/uRNNI_graphs/tree_numbers_%s_taxa.txt' %n, 'output/uRNNI_graphs/uRNNI_edges_%s_taxa.txt' %n)
    else:
    # Compute the uRNNI graph on n taxa
        print("Start generating uRNNI graph on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
        uRNNI = uRNNI_graph_on_n_taxa(n)
        print("Done generating uRNNI graph on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))

    if exists('output/distance_matrices/distance_matrix_%s_taxa.npy' %n):
        D = np.load('output/distance_matrices/distance_matrix_%s_taxa.npy' %n, allow_pickle=False)
        #l = np.load('output/distance_matrices/distance_matrix_index_%s_taxa.npy' %n, allow_pickle=True)
        l = pickle.load(open('output/distance_matrices/distance_matrix_index_%s_taxa.npy' %n, "rb"))
    else:
        print("Start generating distance matrix on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
        S = uRNNI[0].Seidel() #Seidel does not work for more than 7 taxa
        #S.tofile('output/distance_matrices/distance_matrix_%s_taxa' %n)
        np.save('output/distance_matrices/distance_matrix_%s_taxa.npy' %n, S[0], allow_pickle=False)
        #np.save('output/distance_matrices/distance_matrix_index_%s_taxa.npy' %n, S[1], allow_pickle=True)
        pickle.dump(S[1], open('output/distance_matrices/distance_matrix_index_%s_taxa.npy' %n, "wb"), pickle.HIGHEST_PROTOCOL)
        l = S[1]
        D = S[0]


    print("Finding paths & checking their distances on %s at %s" % (time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
    with Pool() as pool:
        equal_distances = p_map(check_tree, uRNNI[1])
        all_correct = len(equal_distances) == sum(equal_distances)
        if (not all_correct):
            return False
    return True

def main():
    n = int(input("For which number of taxa do you want to check FINDPATH?\n"))
    print("Start testing FINDPATH on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
    if check_findpath(n):
        print("Findpath computes correct pairwise distances for all trees on %s taxa" %n)
    else:
        print("Findpath does NOT compute all pairwise distances correctly for %s taxa" %n)
    print("Done testing FINDPATH on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))

if __name__ == "__main__":
    main()
