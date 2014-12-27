CFLAGS = -fopenmp -Wall -O0
CC = gcc-4.9

SOURCES = $(wildcard *.c)
EXES = $(SOURCES:.c=)

all: $(SOURCES) $(EXES)

clean:
	rm -f $(EXES)

%: %.c
	$(CC) $(CFLAGS) $< -o $@
