test_n4:
	mpicc zeta_mpi.c -o z_mpi
	mpirun z_mpi 4
	mpirun z_mpi 8
	mpirun z_mpi 256
	mpirun z_mpi 500000

test_n2:
	mpicc zeta_mpi.c -o z_mpi
	mpirun -n 2 z_mpi 4
	mpirun -n 2 z_mpi 8
	mpirun -n 2 z_mpi 256
	mpirun -n 2 z_mpi 500000

test_n8:
	mpicc zeta_mpi.c -o z_mpi
	mpirun -n 8 z_mpi 4
	mpirun -n 8 z_mpi 8
	mpirun -n 8 z_mpi 256
	mpirun -n 8 z_mpi 500000

test:
	mpicc zeta_mpi.c -o z_mpi
	mpirun z_mpi 2