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
		double a = (double) 1/5;
		double b = (double) 1/239;
		double S1 = 0;
		double S2 = 0;
		double sign = 0;
		for (double i=0;i<=n;i++) {
			sign = pow(-1,i);
			S1 = sign*pow(a,2*i+1)/(2*i+1);
			S2 = sign*pow(b,2*i+1)/(2*i+1);
			S[(int) i] = 4*(4*S1 - S2);
			//printf("%lf\n", S[(int) i]);
		}

	double sum = 0;
	int i;
	int num_proc = 4;
	#pragma omp parallel for reduction(+:sum) num_threads(num_proc) 
	  	for (i=0;i< (int) n;++i) {
			sum += S[i];
		}
	
	double pi = sum;
	printf("pi = %lf\n", pi);
	clock_t diff = clock() - s_time;
	double t_total = (double) (diff) / 1000;
	printf("Absolute error: %.15lf, n = %.0lf, walltime: %.5lfs\n", fabs(PI-pi),n, t_total);

	return 0;
}