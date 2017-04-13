#include "machin_serial.c"


int main(int argc, char const *argv[])
{
	
	if (argc != 2) {
		printf("Expected one argument\n");
		return -1;
	}

	double n = 0;
	n = atof(argv[1]);
	double pi = 0;
	pi = machin0(n);

	printf("pi = %.14lf\n", pi);

	return 0;
}

		

	
