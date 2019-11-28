__author__ = 'Lena Collienne'

#implementation of alg:mdtree
#computes ONE possible tree that could result from T by applying this algorithm
#output does not necessarily have max distance from given tree, e.g. for tree = [set(['1', '2']), set(['1', '2', '4']), set(['3','5']), set(['1', '2', '3', '4', '5'])] it always computes trees with distance 5
from findpath import *
from caterpillar import *
import test_mrcas
import matplotlib.pyplot as plt

def max_dist_tree(T):
    # alternative bottom up procedure to produce caterpillar tree with large distance from given tree T
    n = len(T) + 1
    #First find the order of taxa
    order = dict()
    rank = 1
    has_rank = [False] * n #list showing whether taxon i is already added in order list
    for cluster in T:
        order[rank] = set()
        for element in cluster:
            if has_rank[int(element) - 1] == False:
                order[rank].add(element)
                has_rank[int(element) - 1] = True
        rank += 1
    # Build output tree R
    R = list()
    # The cherry of R shall consist of two taxa a,b on opposite ends of T:
    # a is one of the cherry taxa of T and b is one of the cherry taxa of the lowest cherry on the other side of the root of T
    a = order[1].pop()
    # Find the cluster of the subtree that is child of the root but does not contain a
    S = set()
    for i in range(len(T)-2,-1,-1):
        # print(T[i])
        if a in T[i]:
            S = T[i]
            break
    for i in range(1,len(order)+1):
        if len(order[i]) > 0 and len(S.intersection(order[i])) == 0:
            b = order[i].pop()
            break
    # Build output tree R from bottom to top by reversing the order of the ranks of parents of taxa of T
    R.append(set([a,b]))
    for i in range(len(order),0,-1):
        while (len(order[i]) > 0):
            previous = R.pop()
            R.append(previous)
            R.append(previous.union(set([order[i].pop()])))
    return(R)


def mdtree(T):
    # MDTree as introduced in paper
    # computes a tree with maximum distance from T if T is caterpillar (see paper)
    n = len(T) + 1
    # First find the order of taxa
    order = dict()
    rank = 1
    has_rank = [False] * n # list showing whether taxon i is already added in order list
    for cluster in T:
        order[rank] = set()
        for element in cluster:
            if has_rank[int(element) - 1] == False:
                order[rank].add(element)
                has_rank[int(element) - 1] = True
        rank += 1
    #Build tree
    R = list()
    #Special case: cherry - initial move where both cherry taxa a,b of T appear in last set of R => Now: R = [{a,b}]
    R.append(order[1])
    #Build tree R
    for i in range(2,n):
        while len(order[i]) > 0:
            element = random.choice(tuple(order[i])) #randomly choose element of order[i] to add as next taxon to R
            order[i].remove(element)
            # add new set at first position of R, containing next element in order together with one of the already existing taxa in R (randomly picked out of last set in R)
            x = random.choice(tuple(R[0]))
            for cluster in R:
                if x in cluster:
                    cluster.add(element)
            R.insert(0,set([x,element]))
    return(R)
