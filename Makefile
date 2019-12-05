default:
	gcc -std=c99 -Ofast -fPIC -fopenmp -shared -o libseidel.so seidel.c
