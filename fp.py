__author__ = "Lena Collienne"

from ctypes import *

# Compute RNNI distance between start_tree and end_tree with num_leaves leaves, using FINDPATH implemented in C library (tree.so)
def distance(num_leaves, t1, t2):
    # Convert tree strings into bytes
    start_tree = bytes(t1, 'utf-8')
    end_tree = bytes(t2, 'utf-8')
    lib = cdll.LoadLibrary("./tree.so")
    py_distance = lib.distance
    py_distance.argtypes = [c_int, c_char_p, c_char_p]
    py_distance.restype = c_int
    output = py_distance(num_leaves, start_tree, end_tree)
    return(output)

# Example of execution of distance()
result = distance(5, "[{1,2},{1,2,3},{1,2,3,4},{1,2,3,4,5}]", "[{4,5},{3,4,5},{1,2},{1,2,3,4,5}]")
print('RNNI distance:' , result)