all: libmec.so

CFLAGS=-std=c99 -Wall -pedantic -O3
# LIBS=`pkg-config gsl --libs`
LIBS=-Lmini-gsl -lmini-gsl

libmec.so: mec.c
	gcc -fpic -shared $(CFLAGS) -o $@ $< $(LIBS)
