#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "./TurboPFor/vp4.h"

#include <omp.h>
// TurboPFor likes having a little extra space in it's buffers.
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

void free_matrix(unsigned char **X, int n) {
	for (int i = 0; i < n; i++) {
		free(X[i]);
	}
	free(X);
}

/* Recursive component of Seidel's algorithm.
 * Allocates a new matrix for returning D.
 * Expects the input matrix to be compressed with TurboPFor,
 * and available as both arow & column major version.
 * (for cache-friendly dot-products)
 */
dm seidel_recursive(dm A, int n, int depth) {
	int num_threads = omp_get_max_threads();
	printf("running on %d threads\n", num_threads);
	if (depth > 5) {
		dm D;
		return D;
	}
	size_t matrix_size = (size_t)n * sizeof(int);
	const int verbose = 0; // optionally print entire matrices at each stop for debugging.
	printf("recursion\n");
	// test the base case
	/* int done = 1; */
	int thread_done[num_threads];
	for (int i = 0; i < num_threads; i++) {
		thread_done[i] = 1;
	}

	if (verbose) {
		printf("A:\n");
		print(A, n);
	}
	printf("checking A\n");
	// Allocate separate buffers for each thread.
	uint32_t *row_A[num_threads];
	uint32_t *row_B[num_threads];
	unsigned char *row_B_buf[num_threads];
	uint32_t *col_A[num_threads];
	uint32_t *col_B[num_threads];
	unsigned char *col_B_buf[num_threads];
	for (int i = 0; i < num_threads; i++) {
		row_A[i] = malloc(BUF_SIZE(n)*sizeof(uint32_t));
		row_B[i] = malloc(BUF_SIZE(n)*sizeof(uint32_t));
		row_B_buf[i] = malloc(BUF_SIZE(n)*sizeof(uint32_t));
		col_A[i] = malloc(BUF_SIZE(n)*sizeof(uint32_t));
		col_B[i] = malloc(BUF_SIZE(n)*sizeof(uint32_t));
		col_B_buf[i] = malloc(BUF_SIZE(n)*sizeof(uint32_t));
	}
	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		int thread_id = omp_get_thread_num();
		p4ndec32(A.sa[i], n, row_A[thread_id]);
		for (int j = 0; j < n; j++) {
			if (i != j && row_A[thread_id][j] != 1) {
				thread_done[thread_id] = 0;
			}
		}
	}
	int done = 1;
	for (int i = 0; i < num_threads; i++) {
		if (thread_done[i] == 0)
			done = 0;
	}
	printf("done checking A\n");
	// the array is all ones. We can just return it.
	if (done) {
		printf("floor\n");
		uint32_t *tmp_row_A = malloc(BUF_SIZE(n)*sizeof(uint32_t));
		dm D;
		D.sa = malloc(n*sizeof(unsigned char*));
		for (int i = 0; i < n; i++) {
			size_t row_size = p4ndec32(A.sa[i], n, tmp_row_A);
			D.sa[i] = malloc(row_size);
			memcpy(D.sa[i], A.sa[i], row_size);
		}
		free(tmp_row_A);
		for (int i = 0; i < num_threads; i++) {
			free(row_A[i]);
			free(row_B[i]);
			free(row_B_buf[i]);
			free(col_A[i]);
			free(col_B[i]);
			free(col_B_buf[i]);
		}

		return D;
	}

	printf("calculating B\n");
	dm B;
	uint32_t *degree = malloc(n*sizeof(uint32_t));
	B.sa = malloc(n*sizeof(unsigned char *));
	B.sat = malloc(n*sizeof(unsigned char *));
	// Calculate the row-major matrix B for the next recursion.
	#pragma omp parallel for shared(degree)
	for (int i = 0; i < n; i++) {
		int thread_id = omp_get_thread_num();
		uint32_t deg_i = 0;
		p4ndec32(A.sa[i], n, row_A[thread_id]);
		for (int j = 0; j < n; j++) {
			p4ndec32(A.sat[j], n, col_A[thread_id]);
			deg_i += row_A[thread_id][j];
			if (i == j) {
				row_B[thread_id][j] = 0;
			} else {
				row_B[thread_id][j] = 0;
				if (row_A[thread_id][j] == 1) { // set the only other case where B[i*n+j] would be one while we're here. Avoids the second loop.
					row_B[thread_id][j] = 1;
				}
				for (int k = 0; k < n; k++) {
					// As soon as we've established that Z[i][j] > 0, we can set B[i][j] = 1 and stop.
					// The actual value for Z[i][j] is never used.
					if (row_A[thread_id][k] > 0 && col_A[thread_id][k] > 0) {
						row_B[thread_id][j] = 1;
						break;
					}
				}
			}
		}
		size_t row_b_size = p4nenc32(row_B[thread_id], n, row_B_buf[thread_id]);
		B.sa[i] = malloc(row_b_size);
		memcpy(B.sa[i], row_B_buf[thread_id], row_b_size);

		degree[i] = deg_i;
	}
	printf("done calculating B\n");
	printf("calculating Bt\n");
	// Calculate the column-major matrix B for the next recursion.
	#pragma omp parallel for
	for (int j = 0; j < n; j++) {
		int thread_id = omp_get_thread_num();
		p4ndec32(A.sat[j], n, col_A[thread_id]);
		for (int i = 0; i < n; i++) {
			p4ndec32(A.sa[i], n, row_A[thread_id]);
			if (i == j) {
				col_B[thread_id][i] = 0;
			} else {
				col_B[thread_id][i] = 0;
				if (col_A[thread_id][i] == 1) { // set the only other case where Bt[i*n+j] would be one while we're here.
					col_B[thread_id][i] = 1;
				}
				for (int k = 0; k < n; k++) {
					// As soon as we've established that Z[i][j] > 0, we can set B[i][j] = 1 and stop.
					// The actual value for Z[i][j] is never used.
					if (row_A[thread_id][k] > 0 && col_A[thread_id][k] > 0) {
						col_B[thread_id][i] = 1;
						break;
					}
				}
			}
		}
		size_t col_b_size = p4nenc32(col_B[thread_id], n, col_B_buf[thread_id]);
		B.sat[j] = malloc(col_b_size);
		memcpy(B.sat[j], col_B_buf[thread_id], col_b_size);
	}
	printf("done calculating Bt\n");

	if (verbose) {
		printf("B:\n");
		print(B, n);
	}

	dm T = seidel_recursive(B, n, depth+1);
	if (verbose) {
		printf("T:\n");
		print(T, n);
	}

	printf("calculating D\n");
	free_matrix(B.sa, n);
	free_matrix(B.sat, n);
	dm D;
	D.sa = malloc(n*sizeof(unsigned char*));
	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		int thread_id = omp_get_thread_num();
		// We don't need these buffers anymore, re-use them.
		uint32_t *row_T = row_B[thread_id];
		uint32_t *row_D = col_B[thread_id];
		unsigned char *row_D_buf = row_B_buf[thread_id];
		p4ndec32(T.sa[i], n, row_T);
		for (int j = 0; j < n; j++) {
			p4ndec32(A.sat[j], n, col_A[thread_id]);
			// find vector product of A_{i}, and B_{,j}
			uint32_t prod = 0;
			uint32_t cutoff = row_T[j]*degree[j];
			for (int k = 0; k < n; k++) {
				prod += row_T[k]*col_A[thread_id][k]; // > cutoff exactly when col j of A and row i of T have more than degree[j] overlap.
				if (prod >= cutoff) {
					break;
				}
			}
			row_D[j] = 2*row_T[j];
			if (prod < cutoff)
				row_D[j] -= 1;
		}
		size_t row_d_size = p4nenc32(row_D, n, row_D_buf);
		D.sa[i] = malloc(row_d_size);
		memcpy(D.sa[i], row_D_buf, row_d_size);
	}
	printf("done calculating D\n");

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
	for (int i = 0; i < num_threads; i++) {
		free(row_A[i]);
		free(row_B[i]);
		free(row_B_buf[i]);
		free(col_A[i]);
		free(col_B[i]);
		free(col_B_buf[i]);
	}
	free_matrix(T.sa, n);
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
	size_t matrix_size = (size_t)n * (size_t)n * sizeof(int); //N.B. this will overflow an int.
	unsigned char **sa = malloc(n*sizeof(unsigned char*));
	unsigned int *At = malloc(matrix_size);
	unsigned char **sat = malloc(n*sizeof(unsigned char*));

	unsigned char* comp_A = malloc(BUF_SIZE(n)*sizeof(uint32_t));
	for (int i = 0; i < n; i++) {
		size_t row_a_size = p4nenc32((uint32_t*)&(A[i*n]), n, comp_A);
		uint32_t *buf = malloc((n+127)/128+(n+32)*sizeof(int));
		free(buf);
		sa[i] = malloc(row_a_size);
		memcpy(sa[i], comp_A, row_a_size);
		for (int j = 0; j < n; j++) {
			At[i+j*n] = A[i*n+j];
		}
	}
	for (int i = 0; i < n; i++) {
		size_t col_a_size = p4nenc32(&At[i*n], n, comp_A);
		uint32_t *buf = malloc(BUF_SIZE(n)*sizeof(int));
		sat[i] = malloc(col_a_size);
		memcpy(sat[i], comp_A, col_a_size);
	}
	free(comp_A);
	free(At);
	dm dm_A;
	dm_A.sa = sa;
	dm_A.sat = sat;
	dm D = seidel_recursive(dm_A, n, 0);

	// Convert A into D and free D.
	// This way we clean up everything allocated in C, and don't mess with Python's memory management.
	uint32_t *row_A = malloc(BUF_SIZE(n)*sizeof(uint32_t));
	for (int i = 0; i < n; i++) {
		p4ndec32(D.sa[i], n, row_A);
		for (int j = 0; j < n; j++) {
			A[i*n+j] = row_A[j];
		}
	}
	for (int i = 0; i < n; i++) {
		free(sa[i]);
		free(sat[i]);
	}
	free(sa);
	free(sat);
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
