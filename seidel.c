#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <omp.h>

typedef struct dual_mat {
	short *sa;
	short *sat;
} dm;

void print(dm X, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%d ", X.sat[i+j*n]);
		}
		printf("\n");
	}
}

/* Recursive component of Seidel's algorithm.
 * Allocates a new matrix for returning D.
 * Z doesn't need to be allocated, maybe X does?
 */
void seidel_recursive(short *Dest, dm A, int n, int depth) {
	const int verbose = 0;
	printf("recursion\n");
	// test the base case
	int done = 1;
	if (verbose) {
		printf("A:\n");
		print(A, n);
	}
	for (int i = 0; i < n; i++) {
		#pragma omp parallel for
		for (int j = 0; j < n; j++) {
			if (i != j && A.sa[i*n + j] != 1) {
				done = 0;
			}
		}
	}
	// the array is all ones. We can just return it.
	if (done) {
		printf("floor\n");
		memcpy(Dest, A.sa, n*n*sizeof(short));
		return;
	}
	
	dm B;
	B.sa = malloc(2*n*n*sizeof(short));
	B.sat = malloc(2*n*n*sizeof(short)); //not used
	short *degree = malloc(n*sizeof(short));
	for (int i = 0; i < n; i++) {
		degree[i] = 0;
		#pragma omp parallel for
		for (int j = 0; j < n; j++) {
			degree[i] += A.sa[i*n+j];
			if (i == j) {
				B.sa[i*n+j] = 0;
				B.sat[i+j*n] = 0;
			} else {
				B.sa[i*n+j] = 0;
				B.sat[i+j*n] = 0;
				if (A.sa[i*n+j] == 1) { // set the only other case where B[i*n+j] would be one while we're here.
					B.sa[i*n+j] = 1;
					B.sat[i+j*n] = 1;
				}
				for (int k = 0; k < n; k++) {
					if (A.sa[i*n+k] > 0 && A.sat[k+j*n] > 0) {
						B.sa[i*n+j] = 1;
						B.sat[i+j*n] = 1;
						break;
					}
				}
			}
		}
	}

	if (verbose) {
		printf("B:\n");
		print(B, n);
	}

	seidel_recursive(A.sa, B, n, depth+1);
	dm T;
	T.sa = A.sa;
	free(B.sat); // no longer needed. TODO: also .sa?
	if (verbose) {
		printf("T:\n");
		print(T, n);
	}

	dm D = B;
	for (int i = 0; i < n; i++) {
		#pragma omp parallel for
		for (int j = 0; j < n; j++) {
			// find vector product of A_{i}, and B_{,j}
			short prod = 0;
			short cutoff = T.sa[i*n+j]*degree[j];
			for (int k = 0; k < n; k++) {
				prod += T.sa[i*n+k]*A.sat[k+j*n]; // > cutoff exactly when col j of A and row i of T have more than degree[j] overlap.
				if (prod >= cutoff) {
					break;
				}
			}
			Dest[i*n+j] = 2*T.sa[i*n+j];
			if (prod < cutoff)
				Dest[i*n+j] -= 1;
		}
	}

	if (verbose) {
		printf("degree:\n");
		for (int i = 0; i < n; i++) {
			printf("%d ", degree[i]);
		}
		printf("\n");
	}

	if (verbose) {
		printf("D:\n");
		print(D, n);
	}
	free(T.sa); // is A.sa
	free(degree);

	printf("going up\n");
}


/* A should be an nxn row-major adjacency matrix
 * A is converted in-place into the nxn distance matrix
 * (to avoid passing anything malloc'd back to python.
 * not sure that this is actually a problem, but I don't want to find out)
 */
void seidel(int *A, int n) {
	short *sa = malloc(n*n*sizeof(short));
	short *sat = malloc(n*n*sizeof(short));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			sa[i*n+j] = A[i*n+j];
			sat[i+j*n] = A[i*n+j];
		}
	}
	dm dm_A;
	dm_A.sa = sa;
	dm_A.sat = sat;
	dm D;
	D.sa = malloc(n*n*sizeof(short));
	seidel_recursive(D.sa, dm_A, n, 0);
	//printf("Seidel returned:\n");
	//print(D, n);

	//TODO: convert A into D and free D.
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A[i*n+j] = D.sa[i*n+j];
		}
	}
	free(D.sa);
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
