all: libmec.so test

CFLAGS=-g -std=c99 -Wall -pedantic -O3
LIBS=`pkg-config gsl --libs`
# LIBS=mini-gsl/mini-gsl.o

libmec.so: mec.c
	gcc -fpic -shared $(CFLAGS) -o $@ $< $(LIBS)

test: test.c
	gcc $(CFLAGS) -o $@ $< -L. -lmec -lm
