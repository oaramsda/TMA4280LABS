#include <stdio.h>

#include "machin_serial.c"

#define PI 3.14159265358979323846

int main(int argc, char const *argv[])
{
	FILE *fp;
	fp = fopen("./machin0_residual.txt", "w");

	double pi_k = 0;
	
	for (int k=1; k<=24; ++k) {

		pi_k = machin0(pow(2,k));

		double resi = fabs(PI - pi_k);

		fprintf(fp, "n = %.0lf,\t\t\tresidual = %.14lf\n", pow(2,k), resi);
		

	}

	fclose(fp);
	printf("Done writing to file\n");
}
