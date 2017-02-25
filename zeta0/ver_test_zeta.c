#include <stdio.h>

#include "zeta_serial.c"

#define PI 3.14159265358979323846

int main(int argc, char const *argv[])
{
	FILE *fp;
	fp = fopen("./zeta0_residual.txt", "w");
	
	double pi_k = 0;
	
	for (int k=1; k<=24; ++k) {

		pi_k = zeta0(pow(2,k));

		double resi = PI - pi_k;

		fprintf(fp, "n = %.0lf,\t\t\tresidual = %.15lf, \t\t\tpi=%.15lf\n", pow(2,k), resi, pi_k);
		

	}

	fclose(fp);
	printf("Done writing to file\n");
}