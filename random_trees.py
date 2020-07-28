__author__ = 'Lars Berling'

import rpy2.robjects as robjects
from nwk_parser import *
from summary import *

def r_generate_random_trees(n, k):
    output = []
    rcode = '''
        library(ape)
        t = rmtree(N = {},n = {},rooted=T,tip.label = c(1:{}))
        list = write.tree(phy = t,file = '')
    '''.format(k,n,n)
    list = robjects.r(rcode)
    # list is a string vector ->
    extract = re.compile('\(.*\)\;')
    for t in list :
        cur_tree_str = re.search(extract, t).group()
        cur_tree_str = io.StringIO(cur_tree_str)
        cur_tree = Phylo.read(cur_tree_str, 'newick')
        cur_tree.root.branch_length = float(0)
        rec_nameing(cur_tree.root, cur_tree)
        clust = to_cluster(cur_tree.root.name)
        output.append(string_to_list_of_sets(clust))
    return output
