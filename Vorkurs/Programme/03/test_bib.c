#include "mathe_own.h"
#include <stdio.h>

int   main(){
    double x = 1.35;
    double y = 1.57;
    double z  =  1.;
    printf("Wir testen das eigene Mathemodul und das einbinden:\n");
    printf("Zunächst die sign-Funktion sign(%f) = %d\n", x, sgn_own(x));
    printf("Test des Cosinus cos(%f) = %f\n", y, cos_own(y));
    printf("Test der Exponentialfunktion exp(%f) = %f\n", z, exp_own(z));
    printf("Test des Logarithmus: ln(%f) = %f\n", y, log_own(y));
}