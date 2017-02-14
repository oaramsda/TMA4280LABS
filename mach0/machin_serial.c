#include <stdio.h>
#include <math.h>
#include <stdlib.h>


double machin0(double n)
{
	double S1 = 0;
	double S2 = 0;
	double pi = 0;

	double i = 0;
	double sign = 0;
	double a = (double) 1/5;
	double b = (double) 1/239;

	for (i=0;i<=n;i++) {
		sign = pow(-1,i);
		S1 += sign*pow(a,2*i+1)/(2*i+1);
		S2 += sign*pow(b,2*i+1)/(2*i+1);
	}

	pi = 4*(4*S1 - S2);

	return pi;
}

/*
int main(int argc, char const *argv[])
{
	double n = 0;
	n = atof(argv[1]);
	double pi = 0;
	pi = machin0(n);

	printf("pi = %.14lf\n", pi);

	return 0;
}
*/