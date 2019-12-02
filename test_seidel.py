#!/usr/bin/python3

from uRNNI_graph import *
from os.path import exists
import readline
import math
import cProfile
import re

def check_seidel(n):
    if exists('output/uRNNI_graphs/tree_numbers_%s_taxa.txt' %n):
        uRNNI = read_uRNNI_graph('output/uRNNI_graphs/tree_numbers_%s_taxa.txt' %n, 'output/uRNNI_graphs/uRNNI_edges_%s_taxa.txt' %n)
    else:
        # Compute the uRNNI graph on n taxa
        print("Start generating uRNNI graph on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
        uRNNI = uRNNI_graph_on_n_taxa(n)
        print("Done generating uRNNI graph on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))

    fw_distance = uRNNI[0].Seidel()
    seidel_distance = uRNNI[0].Seidel()

    print(np.array_equal(fw_distance[0], seidel_distance[0]))
    assert(np.array_equal(fw_distance[0], seidel_distance[0]))

def main():
    n = int(input("n: "))
    check_seidel(n)

if __name__ == "__main__":
    #cProfile.run('main()')
    main()
