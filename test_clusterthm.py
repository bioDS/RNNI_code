__author__ = 'Lena Collienne, Kieran Elmes'

#checks whether the cluster thm holds
#only works in reasonable time for up to 7 taxa
#uses symmetry: only checks cluster theorem for trees sharing clusters of the form {1,...,m} for m <= floor(n/2). Result follows for all other clusters through relabelling + symmetry of ranked trees

from uRNNI_graph import *
from os.path import exists
import math

import cProfile
import re


def get_cluster_subset(C, uRNNI):
    #Given a cluster C, this function returns the indices of trees in uRNNI that contain this cluster
    cluster_subset = set()
    n = len(string_to_tree(uRNNI[1][0]))+1 #Take arbitrary tree to define n
    taxa = set() #taxon set {1,..,n}
    for i in range(1,n+1):
        taxa.add(i)
    for i in uRNNI[0].vertices:
        tree = string_to_tree(uRNNI[1][i])
        for node in tree:
            if node == C:
                cluster_subset.add(i)
                break
    return(cluster_subset)


def check_cluster_thm(n):
    # Check if the Cluster Theorem holds by using Seidel on whole uRNNI graph and on the subgraph given by trees containing the same cluster
    if exists('output/uRNNI_graphs/tree_numbers_%s_taxa.txt' %n):
        uRNNI = read_uRNNI_graph('output/uRNNI_graphs/tree_numbers_%s_taxa.txt' %n, 'output/uRNNI_graphs/uRNNI_edges_%s_taxa.txt' %n)
    else:
        # Compute the uRNNI graph on n taxa
        print("Start generating uRNNI graph on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
        uRNNI = uRNNI_graph_on_n_taxa(n)
        print("Done generating uRNNI graph on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))

    C = set()
    C.add(1)
    #only check clusters {1,...,m} for 2 <= m <= floor(n/2); all others follow from symmetry of RNNI
    m = int(math.floor(n/2)) + 1
    uRNNI_distance = uRNNI[0].Seidel()
    uRNNI_adjacency_list = uRNNI[0].adjacency_list()
    pairs_checked = 0
    neighbours_checked = 0
    for i in range(2,m):
        C.add(i)
        print(C)
        # Compute the subset of trees for cluster C in uRNNI graph
        cluster_subset = list(get_cluster_subset(C, uRNNI)) #contains indices of trees in uRNNI containing cluster C
        # Compute the subgraph of uRNNI only containing trees with cluster C + compute all distances within this graph
        print("Start generating cluster graph distance matrix on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
        cluster_graph = uRNNI[0].subgraph(cluster_subset)
        distance_cluster = cluster_graph.Seidel() #Seidel does not work for more than 7 taxa
        print("Done generating cluster graph distance matrix on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
        print("Start comparing cluster %s on %s taxa on %s at %s" % (C, n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))

        for k in cluster_subset:
            for j in cluster_subset:
                pairs_checked += 1
                subset_distance = distance_cluster[0][distance_cluster[1][k]][distance_cluster[1][j]]
                graph_distance = uRNNI_distance[0][uRNNI_distance[1][(k)]][uRNNI_distance[1][(j)]]
                if graph_distance != subset_distance:
                    print("graph distance: ", graph_distance, " != subset_distance: ", subset_distance)
                    print(uRNNI[1][(k)], uRNNI[1][(j)])
                    return False

                for l in uRNNI_adjacency_list[j]:
                    l = int(l)
                    if l not in cluster_subset:
                        neighbours_checked += 1
                        uRNNI_distance[1][l]
                        uRNNI_distance[1][(k)]
                        l_distance = uRNNI_distance[0][uRNNI_distance[1][l]][uRNNI_distance[1][(k)]]
                        j_distance = distance_cluster[0][distance_cluster[1][j]][distance_cluster[1][k]]
                        if l_distance < j_distance:
                            print("l_distance", l_distance)
                            print("j_distance", j_distance)
                            print("l: ", uRNNI[1][l], "is closer than j: ", uRNNI[1][(j)], " to k: ", uRNNI[1][(k)])
                            return False

        print("Done comparing cluster %s on %s taxa on %s at %s" % (C, n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
        print("Checked a total of ", pairs_checked, "pairs, and ", neighbours_checked, "neighbours")
    return True


def main():
    n = int(input("For which number of taxa do you want to check the Cluster Theorem?\n"))
    if check_cluster_thm(n):
        print("Cluster Theorem holds for %s taxa" %n)
    else:
        print("Cluster Theorem does NOT hold for %s taxa" %n)

if __name__ == "__main__":
    main()
