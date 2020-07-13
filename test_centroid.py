__author__ = "Lars Berling"

from summary import *
from tree_generator import *
import timeit

def test_centroid(n, k, x):
    # Function to test if centroid started from all trees in the set does not find a global optimum
    all_trees_on_n_taxa(n)
    print('File Input should be:  all_trees_on_%i_taxa.txt' % n)
    gn = read_trees_from_file()
    for i in range(len(gn)):
        gn[i] = string_to_list_of_sets(gn[i])
    print('The program will check if it will find a case'
          ' where a global optimum is not reached by the centroid algorithm')
    start_time = timeit.default_timer()
    cases = 0
    for i in range(x):
        trees = generate_random_trees(n, k)
        fopt = []
        gopt = global_opt(trees, gn)
        prog = 0
        for s in trees:
            cur_time = timeit.default_timer() - start_time
            prog += 1
            print('Tree %i / %i of set %i / %i after %f minutes' % (prog, k, i + 1, x, cur_time / 60))
            cen = centroid(trees, s)
            fopt.append(min(cen.values()))
        if min(fopt) != gopt:
            cases += 1
        print('---------- Found %i cases so far! ----------' % cases)
    elapsed = timeit.default_timer() - start_time
    print('There were %i cases found where centroid() started from every tree in the set'
          ' does not find a global optimum.\n'
          ' It took %f minutes. ' % (cases, elapsed / 60))
