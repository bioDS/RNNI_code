#include <stdio.h>
#include <string.h>
#include <stdlib.h>


void print(int *X, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%d ", X[i*n+j]);
		}
		printf("\n");
	}
}

int *dot(int *A, int*B, int n) {
	int *R = malloc(n*n*sizeof(int));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			// find vector product of A_{i}, and B_{,j}
			int prod = 0;
			for (int k = 0; k < n; k++) {
				prod += A[i*n+k]*B[k*n+j];
			}
			R[i*n+j] = prod;
		}
	}
	return R;
}

/* Recursive component of Seidel's algorithm.
 * Allocates a new matrix for returning D.
 * Z doesn't need to be allocated, maybe X does?
 */
int *seidel_recursive(int *A, int n, int depth) {
	int verbose = 0;
	printf("recursion\n");
	// test the base case
	int done = 1;
	if (verbose)
		printf("A:\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (verbose)
				printf("%d ", A[i*n + j]);
			if (i != j && A[i*n + j] != 1) {
				done = 0;
			}
		}
		if (verbose)
			printf("\n");
	}
	// the array is all ones. We can just return it.
	if (done) {
		printf("done!\n");
		int *D = malloc(n*n*sizeof(int));
		memcpy(D, A, n*n*sizeof(int));
		return D;
	}
	
	int *Z = malloc(n*n*sizeof(int));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Z[i*n+j] = 0;
			for (int k = 0; k < n; k++) {
				if (A[i*n+k] > 0 && A[k*n+j] > 0) {
					Z[i*n+j] = 1;
					break;
				}
			}
		}
	}
	free(Z);
	Z = dot(A, A, n);

    //B = np.array([
    //    [1 if i != j and (A[i][j] == 1 or Z[i][j] > 0) else 0 for j in range(n)]
    //for i in range(n)])
	// otherwise we need to do stuff
	int *B = malloc(n*n*sizeof(int));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			// find B_{i,j}
			B[i*n+j] = 0;
			if (i != j) {
				if (A[i*n + j] == 1 || Z[i*n+j] > 0) {
					B[i*n+j] = 1;
				}
			}
		}
	}
	free(Z);
	if (verbose) {
		printf("B:\n");
		print(B, n);
	}

	int *T = seidel_recursive(B, n, depth+1);
	if (verbose) {
		printf("T:\n");
		print(T, n);
	}
	free(B);

	int *X = dot(T, A, n);
	if (verbose) {
		printf("X:\n");
		print(X, n);
	}

    //degree = [sum(A[i][j] for j in range(n)) for i in range(n)]
	int *degree = malloc(n*sizeof(int));
	for (int i = 0; i < n; i++) {
		degree[i] = 0;
		for (int j = 0; j < n; j++) {
			degree[i] += A[i*n+j];
		}
	}

	if (verbose) {
		printf("degree:\n");
		for (int i = 0; i < n; i++) {
			printf("%d ", degree[i]);
		}
		printf("\n");
	}

    //D = np.array([
    //    [2 * T[i][j] if X[i][j] >= T[i][j] * degree[j] else 2 * T[i][j] - 1 for j in range(n)]
    //for i in range(n)])
	int *D = malloc(n*n*sizeof(int));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (X[i*n+j] >= T[i*n + j] * degree[j]) {
				D[i*n+j] = 2*T[i*n+j];
			} else {
				D[i*n+j] = 2*T[i*n+j] - 1;
			}
		}
	}
	if (verbose) {
		printf("D:\n");
		print(D, n);
	}
	free(X);
	free(degree);

	return D;
}


/* A should be an nxn row-major adjacency matrix
 * A is converted in-place into the nxn distance matrix
 * (to avoid passing anything malloc'd back to python.
 * not sure that this is actually a problem, but I don't want to find out)
 */
void seidel(int *A, int n) {
	int *D = seidel_recursive(A, n, 0);
	//printf("Seidel returned:\n");
	//print(D, n);

	//TODO: convert A into D and free D.
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A[i*n+j] = D[i*n+j];
		}
	}
	free(D);
}

// A should be an nxn row-major (probably doesn't matter) Adjacency matrix.
int test_function(int *A, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			// row i, position j
			A[i*n+j] = j;
		}
	}

	return n*2;
}
