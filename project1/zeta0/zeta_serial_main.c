#include "zeta_serial.c"

int main(int argc, char const *argv[])
{

	if (argc != 2) {
		printf("Expected one argument\n");
		return -1;
	}

	double n = 0;
	n = atof(argv[1]);
	double pi = 0;

	pi = zeta0(n);
	
	printf("n = %.1lf, pi = %.14lf\n", n, pi);
	return 0;
}