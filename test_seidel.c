#include "seidel.h"
#include <stdlib.h>
#include <stdio.h>

#include "TurboPFor/vp4.h"

int main() {
	int n = 4;
	u_int32_t *A = malloc(n*n*sizeof(int));

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			//if (j == (i + 1) % n || i == (j + 1)%n)
			if (i != j%2)
				A[i*n+j] = 1;
			else
				A[i*n+j] = 0;
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%d ", A[i*n+j]);
		}
		printf("\n");
	}

	int *test = malloc(5*sizeof(int));
	memset(test, 0, 5*sizeof(int));
	for (int i = 0; i < 5; i++)
		test[i] = i;
	for (int i = 0; i < 5; i++) {
		printf("%d\n", test[i]);
	}
	unsigned char *out = malloc((5+127)/128+(5+32)*sizeof(uint32_t));
	size_t sz = p4nenc32(test, 5, out);
	printf("sz: %d\n", sz);
	int *out2 = malloc((5+127)/128+(5+32)*sizeof(uint32_t));
	memset(out2, 0, 5*sizeof(int));
	p4ndec32(out, 5, out2);
	for (int i = 0; i < 5; i++) {
		printf("%d\n", out2[i]);
	}

	printf("&A: %lx\n", A);
	seidel(A, n);
	free(A);
}
