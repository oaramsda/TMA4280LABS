#include <stdio.h>

int main(int argc, char const *argv[])
{
	double n = 10;
	double S = 0;

	printf("Enter n: \n");
	scanf("%lf", &n);
	printf("n = %d\n\n", n);

	double i = 0;
	for (i=1;i<n;i++) {
		S += 1/(i*i);
	}

	printf("S = %.12lf\n", S);
	return 0;
}