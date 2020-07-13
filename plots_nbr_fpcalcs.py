__author__ = "Lars Berling"

from summary import *
import matplotlib.pyplot as plt

def mkdir_p(mypath):
    '''Creates a directory. equivalent to using mkdir -p on the command line'''

    from errno import EEXIST
    from os import makedirs,path

    try:
        makedirs(mypath)
    except OSError as exc: # Python >2.5
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise

def count_fp(n,k):
    # Generates a set of k trees on n taxa and computes the centroid
    # centroid_o() is the function without tracking the already computed paths
    # centroid() is tracking every path that has been computed already
    output = [[],[]]
    trees = generate_random_trees(n,k)
    for i in trees:
        centroid_o(trees,i)
        output[0].append(globals.fp_count)
        centroid(trees,i)
        output[1].append(globals.fp_count)
    return output

def statistic_fpcount(n,k,x):
    # Count the number of find path calculations for k trees on n taxa in x different runs
    # Using the count_fp function to use tracking and non tracking
    output = [[],[]]
    for j in range(x):
        globals.tracking_paths = {} # Resetting the tracking dict after each run --> ?
        cur = count_fp(n, k)
        output[0].extend(cur[0])
        output[1].extend(cur[1])
    return output

def plots_fpcount(n,k,x):
    # Function creates 3 plots comparing the number of findpath calculations for k trees on n taxa and x different sets
    mkdir_p('Plots')
    datapoints = statistic_fpcount(n, k, x)
    fig, box = plt.subplots()
    box.set_title('Boxplot not tracking vs. tracking')
    box.boxplot(datapoints)
    box.set_ylim(0, max(datapoints[1]) + 50)
    box.set_ylabel('#FP calculations')
    plt.savefig('Plots/%i_treesets_size_%i_%i_taxa_boxplot.png' % (x, k, n), bbox_inches='tight')
    # plt.show()
    plt.clf()

    plt.plot(range(len(datapoints[0])), datapoints[0], 'r*', label = 'No tracking')
    plt.plot(range(len(datapoints[1])), datapoints[1], 'b*', label = 'Tracking')
    plt.ylabel('#FP calculations')
    plt.title('Data for %i treesets of size %i on %i taxa' % (x, k, n))
    plt.legend()
    plt.savefig('Plots/%i_treesets_size_%i_%i_taxa.png' % (x, k, n), bbox_inches='tight')
    # plt.show()
    plt.clf()

    plt.plot(range(len(datapoints[1])), datapoints[1], 'b*', label = 'No tracking')
    plt.plot(range(len(datapoints[0])), datapoints[0], 'r*', label = 'Tracking')
    plt.ylabel('#FP calculations')
    plt.ylim(0, max(datapoints[1]) + 50)
    plt.title('Data for %i treesets of size %i on %i taxa, tracking only' % (x, k, n))
    plt.legend()
    plt.savefig('Plots/%i_treesets_size_%i_%i_taxa_tracking.png' % (x, k, n), bbox_inches='tight')
    # plt.show()