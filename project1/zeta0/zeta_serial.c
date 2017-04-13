#include <stdio.h>
#include <math.h>
#include <stdlib.h>


double zeta0(double n)
{
	double S = 0;
	double pi = 0;
	//printf("%lf\n", n);

	double i = 0;
	for (i=1;i<=n;i++) {
		S += 1/(i*i);
	}

	pi = sqrt(6*S);
	return pi;
}