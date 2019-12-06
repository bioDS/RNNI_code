all: seidel seidel_test

seidel:
	gcc -Ofast -std=c99 -fPIC -fopenmp -shared -o libseidel.so seidel.c ./TurboPFor/vp4c.o ./TurboPFor/vp4d.o ./TurboPFor/bitunpack.o ./TurboPFor/bitpack.o ./TurboPFor/bitunpack_sse.o ./TurboPFor/bitutil.o ./TurboPFor/vint.o

seidel_test:
	gcc -g test_seidel.c -o  test_seidel libseidel.so -lm

