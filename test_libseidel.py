#!/usr/bin/env python3

import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer

from uRNNI_graph import *
from os.path import exists

_seidel = ctypes.CDLL("./libseidel.so")
_seidel.test_function.argtypes = (ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int32)
_seidel.seidel.argtypes = (ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int32)
_seidel.seidel_recursive.argtypes = (ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int32, ctypes.c_int32)

def main():
    n = int(input("n: "))
    if exists('output/uRNNI_graphs/tree_numbers_%s_taxa.txt' %n):
        uRNNI = read_uRNNI_graph('output/uRNNI_graphs/tree_numbers_%s_taxa.txt' %n, 'output/uRNNI_graphs/uRNNI_edges_%s_taxa.txt' %n)
    else:
        # Compute the uRNNI graph on n taxa
        print("Start generating uRNNI graph on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))
        uRNNI = uRNNI_graph_on_n_taxa(n)
        print("Done generating uRNNI graph on %s taxa on %s at %s" % (n, time.strftime("%a, %b %d %Y"), time.strftime("%H:%M:%S")))

    AI = uRNNI[0].get_adjacency()
    A = np.ascontiguousarray(AI[0], dtype=np.int32)	
    time1 = time.time()
    _seidel.seidel(A, A.shape[0])
    time2 = time.time()
    print("C Seidel took {:.3f}ms".format((time2 - time1)*1000.0))

    time1 = time.time()
    correct = uRNNI[0].Seidel()[0]
    time2 = time.time()
    print("Py Seidel took {:.3f}ms".format((time2 - time1)*1000.0))
    #correct = uRNNI[0].Floyd_Warshall()[0]
    print("The answer is:")
    print(np.array_equal(correct, A))

if __name__ == "__main__":
    main()
