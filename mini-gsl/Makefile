all: libmini-gsl.so

CFLAGS=-Wall -pedantic -std=c11 -O3

libmini-gsl.so: mini-gsl.c
	gcc -fpic -shared $(CFLAGS) -o $@ $<
