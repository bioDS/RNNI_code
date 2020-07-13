__author__ = "Lars Berling"

from plots_nbr_fpcalcs import *
from global_centroid import *

def plot_size_centroid(n, x):
    # Function that creates boxplots for the size of the centroids (global optima) for x sets on n taxa with different
    # sizes (from 2,...,2n+1)
    datapoints = []
    all_trees_on_n_taxa(n)
    gn = read_trees_from_file('all_trees_on_%i_taxa.txt' % n)
    for i in range(len(gn)):
        gn[i] = string_to_list_of_sets(gn[i])
    for i in range(2, 2 * n + 1):
        datapoints.append(size_of_centroid(n, i, x, gn))
    mkdir_p('Plots')
    fig, box = plt.subplots()
    box.set_title('#Global opts for %i different treesets on %i taxa'% (x,n))
    box.boxplot(datapoints)
    # box.plot([1]*(len(datapoints)+1),'r*')
    box.set_ylabel('#Global Opt Centroids')
    box.set_xlabel('Size of treeset-1')
    plt.savefig('Plots/centroid_size_%i_treesets_%i_taxa_boxplot.png' % (x, n), bbox_inches='tight')
    plt.clf()