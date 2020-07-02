__author__ = "Lars Berling"

from findpath import *
import sys
import globals

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

def update_tracking(path):
    # Update the gloabal dict to track the already computed distances
    # For all trees on the given path, save the distance in the tracking dict
    for i in range(len(path)-1):
        for j in range(i,len(path)):
            if str(expand(path[i])) in globals.tracking_paths:
                # tree i in path has already been visited before
                if str(expand(path[j])) in globals.tracking_paths[str(expand(path[i]))]:
                    # print('ERR: Already computed path gets computed again!')
                    if not globals.tracking_paths[str(expand(path[i]))][str(expand(path[j]))] == j-i:
                        # Check if the distance was already computed and if so check if it is equal
                        print('ERR!: update_tracking() found shortest paths with different lengths!')
                else:
                    globals.tracking_paths[str(expand(path[i]))][str(expand(path[j]))] = j-i
            elif str(path[j]) in globals.tracking_paths:
                # tree j in path has already been visited before
                if str(expand(path[i])) in globals.tracking_paths[str(expand(path[j]))]:
                    # print('ERR: Already computed path gets computed again!')#
                    if not globals.tracking_paths[str(expand(path[j]))][str(expand(path[i]))] == j-i:
                        # Check if the distance was already computed and if so check if it is equal
                        print('ERR!: update_tracking() found shortest paths with different lengths!')
                else:
                    globals.tracking_paths[str(expand(path[j]))][str(expand(path[i]))] = j-i
            else:
                # New entry i.e. neither path[i] nor path[j] were in the tracking dict before
                globals.tracking_paths[str(expand(path[i]))] = {str(expand(path[j])) : j-i}

#####
# update_tracking_m() only puts the start and end node of the path in the tracking dict:
# Improvement similar to saving all trees on the path but not as good --> may be relevant if the memory of the tracking dict gets to big
#####
# def update_tracking_m(path):
#     if str(expand(path[0])) in globals.tracking_paths:
#         if str(expand(path[len(path)-1])) in globals.tracking_paths[str(expand(path[0]))]:
#             if not globals.tracking_paths[str(expand(path[0]))][str(expand(path[len(path)-1]))] == len(path)-1:
#                 print('ERR!')
#         else:
#             globals.tracking_paths[str(expand(path[0]))][str(expand(path[len(path)-1]))] = len(path)-1
#     elif str(expand(path[len(path)-1])) in globals.tracking_paths:
#         if str(expand(path[0])) in globals.tracking_paths[str(expand(path[len(path)-1]))]:
#             if not globals.tracking_paths[str(expand(path[len(path)-1]))][str(expand(path[0]))] == len(path)-1:
#                 print('ERR!')
#         else:
#             globals.tracking_paths[str(expand(path[len(path)-1]))][str(expand(path[0]))] = len(path)-1
#     else:
#         globals.tracking_paths[str(expand(path[0]))] = {str(expand(path[len(path)-1])) : len(path)-1}

def check_tracking(tree1,tree2):
    # Check if the path between the trees tree1 and tree2 has already been computed before, if not return -1
    out = -1
    if str(tree1) in globals.tracking_paths:
        if str(tree2) in globals.tracking_paths[str(tree1)]:
            out = globals.tracking_paths[str(tree1)][str(tree2)]
    elif str(tree2) in globals.tracking_paths:
        if str(tree1) in globals.tracking_paths[str(tree2)]:
            out = globals.tracking_paths[str(tree2)][str(tree1)]
    return out

def get_closer(s, trees, distance):
    # Given a set of trees, a starting vertex and a current distance
    # --> compute trees within one neighbourhood of that vertex that make the distance smaller
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
                # update_tracking_m(cur_p) # Just saving start and end node of the path
                globals.fp_count += 1
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
    globals.fp_count = 0
    # tracking_paths = {}
    startdist = 0
    for i in trees:
        # check tracking_paths 
        check = check_tracking(expand(s),expand(i))
        if check == -1:
            cur_p = findpath(collapse(s), collapse(i))
            path = len(cur_p) - 1
            update_tracking(cur_p)
            # update_tracking_m(cur_p) # Just saving start and end node of the path
            globals.fp_count += 1
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
    # print('#FP Calculations : ', fp_count)
    return (output)


# get_closer_o is a version without keeping track of already computed paths
def get_closer_o(s, trees, distance):
    # Given a set of trees, a starting vertex and a current distance
    # --> compute trees within one neighbourhood of that vertex that make the distance smaller
    Nei = one_neighbourhood(expand(s))
    output = []
    for i in Nei:
        # For each Neibour of start compute trees i that make distance smaller than it is currently
        cur = 0
        for j in trees:
            path = len(findpath(collapse(i), collapse(j)))-1
            globals.fp_count += 1
            cur += path ** 2
            # cur += path
            # maybe only sum of distances -> think if that even makes sense
            # cur += path
        if (cur < distance):
            # Maybe <= but that would be weird, look into that
            output.append([i, cur])
    return (output)

# centroid_o is a version without keeping track of already computed paths
def centroid_o(trees, s):
    # Centroid search starting from s
    # Initialization to run get_closer_function
    globals.fp_count = 0
    startdist = 0
    for i in trees:
        path = len(findpath(collapse(s), collapse(i)))-1
        globals.fp_count += 1
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
    # print('#FP Calculations : ',fp_count)
    return (output)

######### Functions for Checking Centroid() #########
# for a set of trees compute the global minimum sum of squares in the treegraph(whole tree space graph)
def global_opt(trees,treegraph):
    opt = sys.maxsize
    for t in treegraph:
        startdist = 0
        for i in trees:
            check = check_tracking(expand(t), expand(i))
            if check == -1:
                cur_p = findpath(collapse(t), collapse(i))
                path = len(cur_p) - 1
                update_tracking(cur_p)
            else:
                path = check
            startdist += path ** 2
        if startdist < opt:
            opt = startdist
    return(opt)

# for a set of trees compute the global minimum sum of squares in the treegraph(whole tree space graph)
# Version without tracking already computed paths
def slow_global_opt(trees,treegraph):
    opt = sys.maxsize
    for t in treegraph:
        startdist = 0
        for i in trees:
            cur_p = findpath(collapse(t), collapse(i))
            path = len(cur_p) - 1
            startdist += path ** 2
        if startdist < opt:
            opt = startdist
    return(opt)