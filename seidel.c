#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "./TurboPFor/vp4.h"

#include <omp.h>
#define BUF_SIZE(n) ((n+127)/128+(n+32))
typedef struct dual_mat {
	unsigned char **sa;
	unsigned char **sat;
} dm;

void print(dm X, int n) {
	uint32_t *row_A = malloc(BUF_SIZE(n)*sizeof(uint32_t));
	for (int i = 0; i < n; i++) {
		p4ndec32(X.sa[i], n, row_A);
		for (int j = 0; j < n; j++) {
			printf("%d ", row_A[j]);
		}
		printf("\n");
	}
	free(row_A);
}

//unsigned char **allocate_matrix(int n) {
//	unsigned char **X = malloc(n*sizeof(unsigned char*));
//	for (int i = 0; i < n; i++) {
//		X[i] = malloc(n*sizeof(int)); //TODO: this should be smaller, that's the point.
//	}
//	return X;
//}

void free_matrix(unsigned char **X, int n) {
	for (int i = 0; i < n; i++) {
		free(X[i]);
	}
	free(X);
}

/* Recursive component of Seidel's algorithm.
 * Allocates a new matrix for returning D.
 * Z doesn't need to be allocated, maybe X does?
 */
dm seidel_recursive(dm A, int n, int depth) {
	if (depth > 5) {
		dm D;
		return D;
	}
	size_t matrix_size = (size_t)n * sizeof(int);
	const int verbose = 0;
	printf("recursion\n");
	// test the base case
	int done = 1;
	if (verbose) {
		printf("A:\n");
		print(A, n);
	}
	for (int i = 0; i < n; i++) {
		uint32_t *row_A = malloc(BUF_SIZE(n)*sizeof(uint32_t));
		p4ndec32(A.sa[i], n, row_A);
		#pragma omp parallel for
		for (int j = 0; j < n; j++) {
			if (i != j && row_A[j] != 1) {
				done = 0;
			}
		}
		free(row_A);
	}
	// the array is all ones. We can just return it.
	if (done) {
		printf("floor\n");
		uint32_t *row_A = malloc(BUF_SIZE(n)*sizeof(uint32_t));
		dm D;
		D.sa = malloc(n*sizeof(unsigned char*));
		for (int i = 0; i < n; i++) {
			size_t row_size = p4ndec32(A.sa[i], n, row_A);
			D.sa[i] = malloc(row_size);
			memcpy(D.sa[i], A.sa[i], row_size);
		}
		free(row_A);

		return D;
	}

	dm B;
	B.sa = malloc(n*sizeof(unsigned char *));
	B.sat = malloc(n*sizeof(unsigned char *));
	uint32_t *degree = malloc(n*sizeof(uint32_t));
	uint32_t *row_A = malloc(BUF_SIZE(n)*sizeof(uint32_t));
	uint32_t *row_B = malloc(BUF_SIZE(n)*sizeof(uint32_t));
	unsigned char *row_B_buf = malloc(BUF_SIZE(n)*sizeof(uint32_t));
	for (int i = 0; i < n; i++) {
		uint32_t deg_i = 0;
		p4ndec32(A.sa[i], n, row_A);
		#pragma omp parallel for reduction(+: deg_i)
		for (int j = 0; j < n; j++) {
			uint32_t *col_A = malloc(BUF_SIZE(n)*sizeof(uint32_t));
			p4ndec32(A.sat[j], n, col_A);
			deg_i += row_A[j];
			if (i == j) {
				row_B[j] = 0;
			} else {
				row_B[j] = 0;
				if (row_A[j] == 1) { // set the only other case where B[i*n+j] would be one while we're here.
					row_B[j] = 1;
				}
				for (int k = 0; k < n; k++) {
					if (row_A[k] > 0 && col_A[k] > 0) {
						row_B[j] = 1;
						break;
					}
				}
			}
			free(col_A);
		}
		//memcpy(A.sa[i], row_A, n*sizeof(int)); //not changed
		size_t row_b_size = p4nenc32(row_B, n, row_B_buf);
		B.sa[i] = malloc(row_b_size);
		memcpy(B.sa[i], row_B_buf, row_b_size);

		degree[i] = deg_i;
	}
	free(row_A);
	free(row_B);
	free(row_B_buf);
	uint32_t *col_A = malloc(BUF_SIZE(n)*sizeof(uint32_t));
	uint32_t *col_B = malloc(BUF_SIZE(n)*sizeof(uint32_t));
	unsigned char *col_B_buf = malloc(BUF_SIZE(n)*sizeof(uint32_t));
	for (int j = 0; j < n; j++) {
		uint32_t deg_i = 0;
		p4ndec32(A.sat[j], n, col_A);
		#pragma omp parallel for
		for (int i = 0; i < n; i++) {
			uint32_t *row_A = malloc(BUF_SIZE(n)*sizeof(uint32_t));
			p4ndec32(A.sa[i], n, row_A);
			deg_i += row_A[j];
			if (i == j) {
				col_B[i] = 0;
			} else {
				col_B[i] = 0;
				if (col_A[i] == 1) { // set the only other case where B[i*n+j] would be one while we're here.
					col_B[i] = 1;
				}
				for (int k = 0; k < n; k++) {
					if (row_A[k] > 0 && col_A[k] > 0) {
						col_B[i] = 1;
						break;
					}
				}
			}
			free(row_A);
		}
		//memcpy(A.sa[i], row_A, n*sizeof(int)); //not changed
		size_t col_b_size = p4nenc32(col_B, n, col_B_buf);
		B.sat[j] = malloc(col_b_size);
		memcpy(B.sat[j], col_B_buf, col_b_size);
	}
	free(col_B_buf);
	free(col_A);
	free(col_B);

	if (verbose) {
		printf("B:\n");
		print(B, n);
	}

	dm T = seidel_recursive(B, n, depth+1);
	if (verbose) {
		printf("T:\n");
		print(T, n);
	}

	free_matrix(B.sa, n);
	free_matrix(B.sat, n);
	dm D;
	D.sa = malloc(n*sizeof(unsigned char*));
	uint32_t *row_T = malloc(BUF_SIZE(n)*sizeof(uint32_t));
	uint32_t *row_D = malloc(BUF_SIZE(n)*sizeof(uint32_t));
		unsigned char *row_D_buf = malloc(BUF_SIZE(n)*sizeof(int));
	for (int i = 0; i < n; i++) {
		p4ndec32(T.sa[i], n, row_T);
		#pragma omp parallel for
		for (int j = 0; j < n; j++) {
			uint32_t *col_A = malloc(BUF_SIZE(n)*sizeof(uint32_t));
			p4ndec32(A.sat[j], n, col_A);
			// find vector product of A_{i}, and B_{,j}
			uint32_t prod = 0;
			uint32_t cutoff = row_T[j]*degree[j];
			for (int k = 0; k < n; k++) {
				prod += row_T[k]*col_A[k]; // > cutoff exactly when col j of A and row i of T have more than degree[j] overlap.
				if (prod >= cutoff) {
					break;
				}
			}
			row_D[j] = 2*row_T[j];
			if (prod < cutoff)
				row_D[j] -= 1;
			free(col_A);
		}
		size_t row_d_size = p4nenc32(row_D, n, row_D_buf);
		D.sa[i] = malloc(row_d_size);
		memcpy(D.sa[i], row_D_buf, row_d_size);
	}
	free(row_D_buf);
	free(row_T);
	free(row_D);

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
	free_matrix(T.sa, n); // is A.sa
	free(degree);

	printf("going up\n");
	return D;
}


/* A should be an nxn row-major adjacency matrix
 * A is converted in-place into the nxn distance matrix
 * (to avoid passing anything malloc'd back to python.
 * not sure that this is actually a problem, but I don't want to find out)
 */
void seidel(unsigned int *A, int n) {
	//printf("first row:\n");
	//for (int i = 0; i < n; i++) {
	//	printf("%d ", A[i]);
	//}
	//printf("\n");
	size_t matrix_size = (size_t)n * (size_t)n * sizeof(int);
	//unsigned char **sa = allocate_matrix(n);
	unsigned char **sa = malloc(n*sizeof(unsigned char*));
	unsigned int *At = malloc(matrix_size);
	//unsigned char **sat = allocate_matrix(n);
	unsigned char **sat = malloc(n*sizeof(unsigned char*));

	//printf("a\n");
	unsigned char* comp_A = malloc(BUF_SIZE(n)*sizeof(uint32_t));
	for (int i = 0; i < n; i++) {
		size_t row_a_size = p4nenc32((uint32_t*)&(A[i*n]), n, comp_A);
		uint32_t *buf = malloc((n+127)/128+(n+32)*sizeof(int));
		//p4ndec32(comp_A, n, buf);
		//for (int j = 0; j < n; j++) {
		//	printf("%d ", buf[j]);
		//}
		//printf("\n");
		free(buf);
		sa[i] = malloc(row_a_size);
		memcpy(sa[i], comp_A, row_a_size);
		for (int j = 0; j < n; j++) {
			At[i+j*n] = A[i*n+j];
		}
	}
	//printf("b\n");
	for (int i = 0; i < n; i++) {
		size_t col_a_size = p4nenc32(&At[i*n], n, comp_A);
		uint32_t *buf = malloc(BUF_SIZE(n)*sizeof(int));
		//p4ndec32(comp_A, n, buf);
		//for (int j = 0; j < n; j++) {
		//	printf("%d ", buf[j]);
		//}
		//printf("\n");
		sat[i] = malloc(col_a_size);
		memcpy(sat[i], comp_A, col_a_size);
	}
	free(comp_A);
	free(At);
	dm dm_A;
	dm_A.sa = sa;
	dm_A.sat = sat;
	dm D = seidel_recursive(dm_A, n, 0);
	//printf("Seidel returned:\n");
	//print(D, n);

	//convert A into D and free D.
	uint32_t *row_A = malloc(BUF_SIZE(n)*sizeof(uint32_t));
	for (int i = 0; i < n; i++) {
		p4ndec32(D.sa[i], n, row_A);
		for (int j = 0; j < n; j++) {
		//	printf("%d ", row_A[j]);
			A[i*n+j] = row_A[j];
		}
		//printf("\n");
	}
	free(row_A);
	free_matrix(D.sa, n);
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
