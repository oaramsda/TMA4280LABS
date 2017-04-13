hybrid:
	mpicc -fopenmp -o m_hyb mach_hybrid.c
	mpirun m_hyb 500000
all_reduc:
	mpicc zeta_mpi_reduc.c -o z_mpi_reduc
	mpirun -np 8 ./z_mpi_reduc 256
