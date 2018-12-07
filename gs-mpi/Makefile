.PHONY: all

all: gs-subdomain-version

gs-subdomain-version: gs-subdomain-version.c utils.c
	mpicc -O3 -Wall -D TEST_CONV -D PRINT_RESULTS gs-subdomain-version.c utils.c -o gs-subdomain-version