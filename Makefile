CFLAGS=-fopenmp -Wall -O0
CC=gcc
LDLIBS=-lm

SOURCES=$(wildcard *.c)
EXES=$(SOURCES:.c=)

all: $(SOURCES) $(EXES)

clean:
	rm -f $(EXES)
