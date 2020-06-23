__author__ = "Lars Berling"


from caterpillar import *
from findpath import *
import numpy as np

def sort_trees(trees):
    # For a set of trees returned the sorted trees
    output = list()
    for i in trees:
        newtree = []
        for j in i:
            cur = sorted(j)
            newtree.append(cur)
        output.append(newtree)
    return(output)

def expand(tree):
    # Takes a tree in list of lists format, adds leave clusters and returns it as list of sets format
    if(len(tree[0]) == 1):
        return(sort_tree(tree))
    else:        
        tree = add_leaves(tree)
        tree = sort_tree(tree)
        return(tree)

def collapse(tree):
    # Takes a tree and returns list of list format without the leaf cluster
    if(len(tree[0]) == 1):
        output = []
        for i in tree:
            if(len(i) > 1):
                output.append(list(i))
        return(output)
    else:
        return(tree)

################### Algorithm 1 ######################

def get_closer(s,trees,distance):
    # Given a set of trees, a starting vertex and a current distance
    # --> compute trees within one neighbourhood of that vertex that make the distance smaller
    Nei = one_neighbourhood(expand(s))
    output = []
    
    for i in Nei:
        # For each Neibour of start compute trees i that make distance smaller than it is currently
        cur = 0
        for j in trees:
            path = len(findpath(collapse(i),collapse(j)))-1
            cur += path ** 2
            # maybe only sum of distances -> think if that even makes sense
            # cur += path
        if(cur < distance):
            # Maybe <= but that would be weird, look into that
            output.append([i,cur])
    return(output)

def centroid(trees,s):
    # extensive Centroid search starting from s
    # Initialization to run get_closer_function
    startdist = 0
    for i in trees:
        path = len(findpath(collapse(s),collapse(i)))-1
        startdist += path ** 2 
    closer = get_closer(s,trees,startdist)
    output = []
    if not closer:
        output.append([expand(s),startdist])
    else:    
        while bool(closer):
            for i in range(0,len(closer)):
                # for each tree that is closer
                cur = get_closer(closer[i][0],trees,closer[i][1])    
                if bool(cur):
                    # there are trees that give a smaller distance = cur
                    # we need to add each tree to close, i.e. check if there is anything in their neighbourhood
                    # which gives better distance
                    for j in cur:
                        closer.append(j)
                else:
                    # No tree in one neighbourhood of i is getting any better distance
                    # --> i is a local optimum and gets put into the output and deleted from things to work on
                    output.append(closer[i])
                closer.pop(i)
    return(output)

