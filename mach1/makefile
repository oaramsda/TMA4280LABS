test_n4:
	mpicc mach_mpi.c -o m_mpi
	mpirun m_mpi 4
	mpirun m_mpi 8
	mpirun m_mpi 256
	mpirun m_mpi 500000

test_n2:
	mpicc mach_mpi.c -o m_mpi
	mpirun -np 2 m_mpi 4
	mpirun -np 2 m_mpi 8
	mpirun -np 2 m_mpi 256
	mpirun -np 2 m_mpi 500000

test_n8:
	mpicc mach_mpi.c -o m_mpi
	mpirun -np 8 m_mpi 4
	mpirun -np 8 m_mpi 8
	mpirun -np 8 m_mpi 256
	mpirun -np 8 m_mpi 500000

test:
	mpicc mach_mpi.c -o m_mpi
	mpirun m_mpi 2