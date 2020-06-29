__author__ = "Lars Berling"

# setwd("~/Desktop/CodingMA/RNNI_code")

from findpath import *

global fp_count
global tracking_paths

def sort_trees(trees):
    # For a set of trees returned the sorted trees
    output = list()
    for i in trees:
        newtree = []
        for j in i:
            cur = sorted(j)
            newtree.append(cur)
        output.append(newtree)
    return (output)


def expand(tree):
    # Takes a tree in list of lists format, adds leave clusters and returns it in list of sets format
    if (len(tree[0]) == 1):
        return (sort_tree(tree))
    else:
        tree = add_leaves(tree)
        tree = sort_tree(tree)
        return (tree)


def collapse(tree):
    # Takes a tree and returns list of list format without the leaf cluster
    if (len(tree[0]) == 1):
        output = []
        for i in tree:
            if (len(i) > 1):
                output.append(list(i))
        return (output)
    else:
        return (tree)
################### Alex Algorithm ######################
# Go from 1 tree and consider a tree within the 1-neighbour hood that minimized the sum of distances t oll other trees given
# Do this until no more can be done
def update_tracking(path):
    global tracking_paths
    for i in range(len(path)-1):
        for j in range(i,len(path)):
            if str(expand(path[i])) in tracking_paths:
                # tree i in path has already been visited before
                if str(expand(path[j])) in tracking_paths[str(expand(path[i]))]:
                    # print('ERR: Already computed path gets computed again!')
                    if not tracking_paths[str(expand(path[i]))][str(expand(path[j]))] == j-i:
                        print('ERR!')
                else:
                    tracking_paths[str(expand(path[i]))][str(expand(path[j]))] = j-i
            elif str(path[j]) in tracking_paths:
                # tree j in path has already been visited before
                if str(expand(path[i])) in tracking_paths[str(expand(path[j]))]:
                    # print('ERR: Already computed path gets computed again!')#
                    if not tracking_paths[str(expand(path[j]))][str(expand(path[i]))] == j-i:
                        print('ERR!')
                else:
                    tracking_paths[str(expand(path[j]))][str(expand(path[i]))] = j-i
            else:
                tracking_paths[str(expand(path[i]))] = {str(expand(path[j])) : j-i}

def check_tracking(tree1,tree2):
    global tracking_paths
    out = -1
    if str(tree1) in tracking_paths:
        if str(tree2) in tracking_paths[str(tree1)]:
            out = tracking_paths[str(tree1)][str(tree2)]
    elif str(tree2) in tracking_paths:
        if str(tree1) in tracking_paths[str(tree2)]:
            out = tracking_paths[str(tree2)][str(tree1)]
    return out


def get_closer(s, trees, distance):
    # Given a set of trees, a starting vertex and a current distance
    # --> compute trees within one neighbourhood of that vertex that make the distance smaller
    global fp_count
    global tracking_paths
    Nei = one_neighbourhood(expand(s))
    output = []
    for i in Nei:
        # For each Neibour of start compute trees i that make distance smaller than it is currently
        cur = 0
        for j in trees:
            # check tracking_paths
            # check expand(i) and expand(j) is already in tracking_paths
            check = check_tracking(expand(i),expand(j))
            if check == -1:
                cur_p = findpath(collapse(i), collapse(j))
                path = len(cur_p) - 1
                update_tracking(cur_p)
                fp_count += 1
            else:
                path = check
            cur += path ** 2
            # cur += path
            # maybe only sum of distances -> think if that even makes sense
            # cur += path
        if (cur < distance):
            # Maybe <= but that would be weird, look into that
            output.append([i, cur])
    return (output)


def centroid(trees, s):
    # Centroid search starting from s
    # Initialization to run get_closer_function
    global fp_count
    global tracking_paths
    fp_count = 0
    tracking_paths = {}
    startdist = 0
    for i in trees:
        # check tracking_paths 
        check = check_tracking(expand(s),expand(i))
        if check == -1:
            cur_p = findpath(collapse(s), collapse(i))
            path = len(cur_p) - 1
            update_tracking(cur_p)
            fp_count += 1
        else:
            path = check
        startdist += path ** 2
    closer = get_closer(s, trees, startdist)
    output = {}
    if not closer:
        output[str(expand(s))] = startdist
    else:
        while bool(closer):
            update = []
            for i in range(len(closer)):
                # for each tree that is closer
                cur = get_closer(closer[i][0], trees, closer[i][1])
                if bool(cur):
                    # there are trees that give a smaller distance = cur
                    # we need to add each tree to closer, i.e. check if there is anything in their neighbourhood
                    # which gives better distance
                    for j in cur:
                        update.append(j)
                else:
                    # No tree in one neighbourhood of i is getting any better distance
                    # --> i is a local optimum and gets put into the output and deleted from things to work on
                    if not str(closer[i][0]) in output:
                        output[str(closer[i][0])] = closer[i][1]
                    # else:
                    # This would be the case if there are multiple paths that are considered to get
                    # to this optimal node
                    # print('Err: Found an interesting example?')
                #closer.pop(i)
            closer = update
    print('#FP Calculations : ', fp_count)
    return (output)


def get_closer_o(s, trees, distance):
    # Given a set of trees, a starting vertex and a current distance
    # --> compute trees within one neighbourhood of that vertex that make the distance smaller
    global fp_count
    Nei = one_neighbourhood(expand(s))
    output = []
    for i in Nei:
        # For each Neibour of start compute trees i that make distance smaller than it is currently
        cur = 0
        for j in trees:
            path = len(findpath(collapse(i), collapse(j)))-1
            fp_count += 1
            cur += path ** 2
            # cur += path
            # maybe only sum of distances -> think if that even makes sense
            # cur += path
        if (cur < distance):
            # Maybe <= but that would be weird, look into that
            output.append([i, cur])
    return (output)


def centroid_o(trees, s):
    # Centroid search starting from s
    # Initialization to run get_closer_function
    global fp_count
    fp_count = 0
    startdist = 0
    for i in trees:
        path = len(findpath(collapse(s), collapse(i)))-1
        fp_count += 1
        startdist += path ** 2
    closer = get_closer_o(s, trees, startdist)
    output = {}
    if not closer:
        output[str(expand(s))] = startdist
    else:
        while bool(closer):
            update = []
            for i in range(len(closer)):
                # for each tree that is closer
                cur = get_closer_o(closer[i][0], trees, closer[i][1])
                if bool(cur):
                    # there are trees that give a smaller distance = cur
                    # we need to add each tree to closer, i.e. check if there is anything in their neighbourhood
                    # which gives better distance
                    for j in cur:
                        update.append(j)
                else:
                    # No tree in one neighbourhood of i is getting any better distance
                    # --> i is a local optimum and gets put into the output and deleted from things to work on
                    if not str(closer[i][0]) in output:
                        output[str(closer[i][0])] = closer[i][1]
                    # else:
                    # This would be the case if there are multiple paths that are considered to get
                    # to this optimal node
                    # print('Err: Found an interesting example?')
                #closer.pop(i)
            closer = update
    print('#FP Calculations : ',fp_count)
    return (output)
