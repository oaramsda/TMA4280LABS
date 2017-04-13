#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define PI 3.14159265358979323846

int main(int argc, char const *argv[])
{
	clock_t s_time = clock();

	MPI_Init(NULL, NULL);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	
	//printf("%lf\n", world_size % 2);
	double l = log(world_size)/log(2);

	if (floor(l) != l) {
		if (world_rank == 0) {
			printf("Number of processes is not a power of 2. Aborting\n");
		}
		MPI_Finalize();
		return 0;
	}

	if (argc != 2) {
		printf("Expected one argument\n");
		MPI_Finalize();
		return -1;
	}
	double n = atof(argv[1]);
	int k = ((int) n) / world_size;
	double *S;
	double *s;
	double par_sum = 0;
	s = (double*)calloc(k, sizeof(double));

	if (world_rank == 0) {
		S = (double*)calloc(n, sizeof(double));
		for (double i=1;i<=n;i++) {
			S[(int) i] = 1/(i*i);
		}
	}
	
	MPI_Scatter(S, k, MPI_DOUBLE, s, k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	for (int i=0;i<k;++i){
		par_sum += s[i];
	}
	if (world_rank == 0) {
		for (int i=k*world_size;i<n;++i){
			par_sum += S[i];
		}
	}

	double sum = 0;
	MPI_Allreduce(&par_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	double pi = sqrt(6*sum);
	printf("pi = %lf\n", pi);
	clock_t diff = clock() - s_time;
	double t_total = (double) (diff) / 1000;
	printf("Absolute error: %.15lf, n = %.0lf, walltime: %.3lfs, process: %d\n", fabs(PI-pi),n, t_total, world_rank);

	
	MPI_Finalize();


	return 0;
}