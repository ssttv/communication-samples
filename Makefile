CFLAGS = -fopenmp -Wall -O0
CC = gcc
LIBS=-lm

SOURCES = $(wildcard *.c)
EXES = $(SOURCES:.c=)

all: $(SOURCES) $(EXES)

clean:
	rm -f $(EXES)

%: %.c
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)
