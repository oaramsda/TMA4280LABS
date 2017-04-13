CC=gcc

run:
	$(CC) -fopenmp -o m_openmp mach_openmp.c
	./m_openmp 4

test:
	$(CC) -fopenmp -o m_openmp mach_openmp.c
	./m_openmp 4
	./m_openmp 8
	./m_openmp 256
	./m_openmp 500000