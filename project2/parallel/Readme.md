How to compile and run:

- CC=mpicc FC=mpif90 cmake ./src -DCMAKE_BUILD_TYPE=Release

- make

- mpirun -np 4 ./poisson 128 2
  - where 4 is number of processes, 128 is problem/grid size, 2 number of threads per process.
