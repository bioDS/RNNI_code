__author__ = 'Lena Collienne'

#checks whether FINDPATH computes correct distances
#works on small trees with up to 6 taxa

from findpath import *
from os.path import exists
import re


def check_findpath(n):
    # Check if FINDPATH computes exact distances by comparing these distances with exact once (Floyd Warshall on RNNI graph)
    if exists('output/uRNNI_graphs/tree_numbers_%s_taxa.txt' %n):
        uRNNI = read_uRNNI_graph('output/uRNNI_graphs/tree_numbers_%s_taxa.txt' %n, 'output/uRNNI_graphs/uRNNI_edges_%s_taxa.txt' %n)
    else:
    # Compute the uRNNI graph on n taxa
        print("Start generating uRNNI graph on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
        uRNNI = uRNNI_graph_on_n_taxa(n)
        print("Done generating uRNNI graph on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))

    FW = uRNNI[0].Floyd_Warshall() #FW does not work for more than 6 taxa
    D = FW[0]
    l = FW[1]

    count = 0
    progress = 0
    max_count = len(uRNNI[1])**2
    for tree1 in uRNNI[1]:
        for tree2 in uRNNI[1]:
            count += 1
            if (len(findpath(string_to_tree(str(tree1)),string_to_tree(str(tree2))))-1 != D[l[uRNNI[1].index(tree1)]][l[uRNNI[1].index(tree2)]]):
                print(tree1, tree2)
                return False
            # Print progress
            if((count/max_count) > progress):
                    print ('Progress: ', str(100*progress), '% at', time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S"))
                    progress += 0.01
    print(count)
    return True


def main():
    n = int(input("For which number of taxa do you want to check FINDPATH?\n"))
    print("Start testing FINDPATH on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
    if check_findpath(n):
        print("Findpath computes correct pairwise distances for all trees on %s taxa" %n)
    else:
        print("Findpath does NOT ocmpute all pairwise distances correctly for %s taxa" %n)
    print("Done testing FINDPATH on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))

if __name__ == "__main__":
    main()
