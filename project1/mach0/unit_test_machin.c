#include "../minunit.h"
#include "machin_serial.c"


int tests_run = 0;
double pi3 = 0;
double pi_3 = 0;
char message[50];

static char *test_machin0()
{	
	double a = (double) 1/5;
	double b = (double) 1/239;
	double S1 = a - pow(a,3)/(3) + pow(a,5)/(5) - pow(a,7)/(7);
	double S2 = b - pow(b,3)/(3) + pow(b,5)/(5) - pow(b,7)/(7);
	pi_3 = 4*(4*S1-S2);
	pi3 = machin0(3.0);
	sprintf(message,"Error, EXPECTED: machin0(3) = %.6lf... WAS: machin0(3) = %lf", pi_3, pi3);
	mu_assert(message, pi3 == pi_3);
	return 0;
}

static char *all_tests() 
{
	mu_run_test(test_machin0);
	return 0;
}

int main(int argc, char const *argv[])
{
	char *result = all_tests();
    if (result != 0) {
         printf("%s\n", result);
    }
    else {
        printf("ALL TESTS PASSED\n");
    }
    printf("Tests run: %d\n", tests_run);
 
    return result != 0;
}