__author__ = "Lars Berling"
from summary import *

def sos_dict(trees,treegraph):
    # returns a dict with the sum of squares for a fixed set of trees and all trees in treespace
    output = {}
    for t in treegraph:
        startdist = 0
        for i in trees:
            cur_p = findpath(collapse(t), collapse(i))
            path = len(cur_p) - 1
            startdist += path ** 2
        # print('The cur sum of squares is %i'%startdist)
        output[str(collapse(t))] = startdist
    return(output)

def get_min_keys(d):
    # returns a list of the keys corresponding to the minimal value in a dict
    mv = min(d.values())
    res = list(filter(lambda x: d[x] == mv, d))
    return(res)

def is_set_connected(trees):
    # Returns true if the set of trees is connected, otherwise false
    if len(trees) == 1:
        return True
    # Initialization
    treeset = set()
    for i in trees:
        treeset.add(str(expand(i)))
    s = trees[0]
    comp_set = set([str(expand(trees[0]))])
    sn = one_neighbourhood(expand(s))
    for i in range(len(sn)):
        sn[i] = str(sn[i])
    sn = set(sn)
    working_list = sn.intersection(set(treeset))
    # End of initialization
    if not working_list:
        # Starting tree is already isolated from the rest of the trees in trees
        return False
    # Starting while loop
    while working_list:
        cur = working_list.pop()
        curnei = one_neighbourhood(expand(string_to_list_of_sets(cur)))
        for i in range(len(curnei)):
            curnei[i] = str(curnei[i])
        curnei = set(curnei)
        rel = curnei.intersection(set(treeset))
        rel = rel.difference(comp_set)
        working_list.update(rel)
        comp_set.add(cur)
    if len(trees) == len(comp_set):
        return True
    return False

def is_centroid_connected(n, k, x):
    # Function that computes all global optimal centroids for n taxa, k trees and x random sets
    # and tests if they are connected or not -> retruns a list of bools
    # If there is only one centroid we enter 'Unique Centroid' to the output --> object to change
    all_trees_on_n_taxa(n)
    print('File Input should be all_trees_on_%i_taxa.txt'%n)
    gn = read_trees_from_file()
    for i in range(len(gn)):
        gn[i] = string_to_list_of_sets(gn[i])
    out = [False]*x
    for i in range(x):
        trees = generate_random_trees(n,k)
        sos = sos_dict(trees,gn)
        mk = get_min_keys(sos)
        mcen = []
        if len(mk) > 1:
            for j in mk:
                mcen.append(string_to_list_of_sets(j))
            if is_set_connected(mcen):
                out[i] = True
            # else : would give us an unconnected centroid
        else:
            # Case that there is only one global opt/centroid --> seems to be very unlikely with 6 taxa
            out[i] = 'Unique Centroid'
    return(out)

def size_of_centroid(n,k,x):
    # Function that computes the global optimal solution for n taxa, k trees and x random sample sets and returns a list
    # with the size of the centroid / how many trees minimize the sum of squared distances globally
    out = []
    all_trees_on_n_taxa(n)
    gn = read_trees_from_file('all_trees_on_%i_taxa.txt'%n)
    for i in range(len(gn)):
        gn[i] = string_to_list_of_sets(gn[i])
    for i in range(x):
        trees = generate_random_trees(n, k)
        sos = sos_dict(trees, gn)
        mk = get_min_keys(sos)
        out.append(len(mk))
    return(out)