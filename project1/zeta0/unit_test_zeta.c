#include <stdio.h>
#include <math.h>

#include "../minunit.h"
#include "zeta_serial.c"

int tests_run = 0;
double pi3 = 0;
double pi_3 = 0;
char message[50];

static char *test_zeta0()
{
	pi_3 = sqrt(6*(1 + (double) 1/(2*2) + (double) 1/(3*3)));
	pi3 = zeta0(3.0);
	sprintf(message,"Error, EXPECTED: zeta0(3) = %.6lf... WAS: zeta0(3) = %lf", pi_3, pi3);
	mu_assert(message, pi3 == pi_3);
	return 0;
}

static char *all_tests() 
{
	mu_run_test(test_zeta0);
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