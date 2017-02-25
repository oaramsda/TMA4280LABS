test_n4:
	mpicc zeta_mpi.c -o z_mpi
	mpirun z_mpi 2
	mpirun z_mpi 4
	mpirun z_mpi 8
	mpirun z_mpi 16
	mpirun z_mpi 32
	mpirun z_mpi 64
	mpirun z_mpi 128
	mpirun z_mpi 256
	mpirun z_mpi 32000
	mpirun z_mpi 124000
	mpirun z_mpi 496000

test_n2:
	mpicc zeta_mpi.c -o z_mpi
	mpirun -n 2 z_mpi 2
	mpirun -n 2 z_mpi 4
	mpirun -n 2 z_mpi 8
	mpirun -n 2 z_mpi 16
	mpirun -n 2 z_mpi 32
	mpirun -n 2 z_mpi 64
	mpirun -n 2 z_mpi 128
	mpirun -n 2 z_mpi 256
	mpirun -n 2 z_mpi 32000
	mpirun -n 2 z_mpi 124000
	mpirun -n 2 z_mpi 496000

test_n8:
	mpicc zeta_mpi.c -o z_mpi
	mpirun -n 8 z_mpi 2
	mpirun -n 8 z_mpi 4
	mpirun -n 8 z_mpi 8
	mpirun -n 8 z_mpi 16
	mpirun -n 8 z_mpi 32
	mpirun -n 8 z_mpi 64
	mpirun -n 8 z_mpi 128
	mpirun -n 8 z_mpi 256
	mpirun -n 8 z_mpi 32000
	mpirun -n 8 z_mpi 124000
	mpirun -n 8 z_mpi 496000

test:
	mpicc zeta_mpi.c -o z_mpi
	mpirun z_mpi 2