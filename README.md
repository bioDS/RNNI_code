# Geometry of Ranked Nearest Neighbour Interchange Space of Phylogenetic Trees - Code Repository

This repository contains the code used to check the FINDPATH algorithm, and the Cluster Theorem from the paper [Geometry of Ranked Nearest Neighbour Interchange Space of Phylogenetic Trees](https://doi.org/10.1101/2019.12.19.883603) by Lena Collienne, Kieran Elmes, Mareike Fischer, David Bryant, and Alex Gavryushkin.


## Compilation
Requires gcc with OpenMP support.

`git submodule init`

`git submodule update`

`make`


## Files
File			|	Description
---			|	---
findpath.py		|	FINDPATH algorithm
test_findpath.py	|	Checks that FINDPATH finds the correct distances for every pair of trees on up to n taxa.
				Relies on POSIX fork behaviour for multiprocessing, probably doesn't work on Windows.
				Distance matrices calculated using Seidel's Algorithm are stored in ./output/distance_matrices
test_clusterthm.py	|	Checks that the Cluster Theorem is true by comparing (for all pairs of trees) the true distance with the distance in the subgraph that maintains the cluster.
			 	If the distance of the path found by FINDPATH differs from the true distance between any pair of trees, we print the offending trees and fail.
				(Distance matrices are not stored/loaded here).
				For every pair of trees S and T, where the distance from S to T is k, we check that the distance from every neighbour S' of S to T is only k-1 if S' maintains the cluster.
				If this condition is not true for pair of trees we print the offending tree and fail.
graph.py		|	Contains a collection of useful algorithms for the RNNIu graph (Seidel, Floyd-Warshall, breadth-first search shortest path, Djikstra, Subgraph defined by a set of vertices).
seidel.c		|	C implementation of Seidel's algorithm.
				Calling seidel() in graph.py uses this version.
				Relies on OpenMP for parallelization.
				Matrices are compressed using TurboPFor, which is included as a git submodule.
				Runs on a 56700x56700 adjacency matrix in ~30 hours, and using ~20GB of memory (on 4x Opteron 6276).
				To actually use this (and by extension, test_findpath or test_clusterthm), first run `git submodule init; git submodule update; make`
