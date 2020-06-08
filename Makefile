CC? = gcc

all: TurboPFor/vp4d.o seidel 

seidel:
	${CC} -Ofast -std=c99 -fPIC -fopenmp -shared -o libseidel.so seidel.c ./TurboPFor/vp4c.o ./TurboPFor/vp4d.o ./TurboPFor/bitunpack.o ./TurboPFor/bitpack.o ./TurboPFor/bitunpack_sse.o ./TurboPFor/bitutil.o ./TurboPFor/vint.o

seidel_test:
	${CC} -g test_seidel.c -o  test_seidel libseidel.so -lm

TurboPFor/vp4d.o:
	CFLAGS=-fPIC make -C TurboPFor vp4d.o vp4c.o bitunpack.o bitpack.o bitutil.o vint.o bitunpack_sse.o

clean:
	find . -type f -name libseidel.o -delete -or -name test_seidel -delete
	make -C TurboPFor clean
