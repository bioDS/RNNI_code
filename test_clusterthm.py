__author__ = 'Lena Collienne'

#checks whether the cluster thm holds
#only works in reasonable time for up to 6 taxa
#uses symmetry: only checks cluster theorem for trees sharing clusters of the form {1,...,m} for m <= floor(n/2). Result follows for all other clusters through relabelling + symmetry of ranked trees

from uRNNI_graph import *
from os.path import exists
import math


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
    # Check if the Cluster Theorem holds by using Dijkstra on whole uRNNI graph and Floyd Warshall on subgraph given by trees containing the same cluster
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
    for i in range(2,m):
        C.add(i)
        print(C)
        # Compute the subset of trees for cluster C in uRNNI graph
        cluster_subset = list(get_cluster_subset(C, uRNNI)) #contains indices of trees in uRNNI containing cluster C
        # Compute the subgraph of uRNNI only containing trees with cluster C + compute all distances within this graph
        print("Start generating cluster graph distance matrix on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
        cluster_graph = uRNNI[0].subgraph(cluster_subset)
        distance_cluster = cluster_graph.Floyd_Warshall() #FW does not work for more than 6 taxa
        print("Done generating cluster graph distance matrix on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
        print("Start comparing cluster %s on %s taxa on %s at %s" % (C, n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
        uRNNI_distance = np.array(np.zeros((len(cluster_subset),len(cluster_subset)))) #exact RNNI distances between trees with shared cluster
        Q = list()
        for k in cluster_subset:
            for j in cluster_subset:
                Q.append([k,j])
        count = 0
        max_count = len(Q)
        h = 0
        while len(Q)>0:
            if(1-len(Q)/max_count > h):
                print ('Progress cluster ', C, ': ', 100*h, '% at', time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S"))
                h += 0.01
            pair = Q.pop()
            l = pair[0]
            D = uRNNI[0].Dijkstra(l)[0] #dict: tree numbers as keys, distance to tree i as values
            for k in D:
                if [k,l] in Q:
                    Q.remove([k,l])
                    Q.remove([l,k])
                    uRNNI_distance[cluster_subset.index(l)][cluster_subset.index(k)] = D[k]
                    uRNNI_distance[cluster_subset.index(k)][cluster_subset.index(l)] = D[k]

        for k in cluster_subset:
            for j in cluster_subset:
                if uRNNI_distance[cluster_subset.index(k)][cluster_subset.index(j)] != distance_cluster[0][distance_cluster[1][k]][distance_cluster[1][j]]:
                    print(uRNNI[1][cluster_subset.index(k)], uRNNI[1][cluster_subset.index(j)])
                    return False
        print("Done comparing cluster %s on %s taxa on %s at %s" % (C, n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
    return True


def main():
    n = int(input("For which number of taxa do you want to check the Cluster Theorem?\n"))
    if check_cluster_thm(n):
        print("Cluster Theorem holds for %s taxa" %n)
    else:
        print("Cluster Theorem does NOT hold for %s taxa" %n)

if __name__ == "__main__":
    main()
