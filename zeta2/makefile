CC=gcc

run:
	$(CC) -fopenmp -o z_openmp zeta_openmp.c
	./z_openmp 5000000

test:
	$(CC) -fopenmp -o z_openmp zeta_openmp.c
	./z_openmp 4
	./z_openmp 8
	./z_openmp 256
	./z_openmp 500000