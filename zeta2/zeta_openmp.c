#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define PI 3.14159265358979323846

int main(int argc, char const *argv[])
{
	clock_t s_time = clock();

	double n = atof(argv[1]);

	double *S;

	S = (double*)calloc(n, sizeof(double));
	for (double i=1;i<=n;i++) {
		S[(int) i] = 1/(i*i);
		//printf("%lf\n", S[(int) i]);
	}

	double sum = 0;
	int i;
	#pragma omp parallel for reduction (+:sum)
	  	for (i=0;i< (int) n;++i) {
			sum += S[i];
		}
	
	double pi = sqrt(6*sum);
	printf("pi = %lf\n", pi);
	clock_t diff = clock() - s_time;
	double t_total = (double) (diff) / CLOCKS_PER_SEC;
	printf("Absolute error: %.10lf, n = %.0lf, walltime: %.5lfs\n", fabs(PI-pi),n, t_total);

	return 0;
}