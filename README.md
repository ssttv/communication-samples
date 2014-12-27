# Communication in parallel applications based on shared memory: sample code


This project contains a set of small benchmarks to show various different communication patterns in parallel applications based on the shared memory paradigm (more specifically, the OpenMP standard).


## Benchmarks

The following benchmarks are provided:

### init
Simple initialization code with an initialization communication pattern.


### reduce
Simple reduction code with a reduction communication pattern.


### matrix
Matrix multiplication with an all-to-all communication pattern.


### gauss-seidel
Gauss-Seidel method with a nearest-neighbor communication pattern.


### gauss-blur
Gaussian blur with a matrix communication pattern.


### prodcons
Producer-consumer benchmark with a pipeline communication pattern.



## Compilation

Compile the benchmarks with

    make

## Execution

Benchmarks inputs have been chosen such that the patterns are best visible with 8 threads.
Therefore, executing with 8 threads is recommended, even on platforms that support executing more or less than 8 threads at the same time.


    export OMP_NUM_THREADS=8
    ./init

