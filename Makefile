all: libmec.so

CFLAGS=-std=c99 -Wall -pedantic -O3

libmec.so: mec.c
	gcc -fpic -shared $(CFLAGS) -o $@ $< `pkg-config gsl --libs`
