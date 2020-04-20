#include "mathe_own.h"
#include <stdio.h>

double root_to(double* a){
    double x = (*a);
    x = wurzel_own(x);
    (*a) = x;
    return x;
}

double square_to(double* a){
    double x = (*a);
    x += x;
    (*a) = x;
    return x;
}

void square_equ(double a, double b, double c){
    double p = b/a;
    double q = c/a;
    if((p*p/4 - q) < 0){
        printf("Es gibt keine Lösungen.\n");
        return;
    }
    double x1 = -p/2 + sqrt( p*p/4. - q );
    double x2 = -p/2 - wurzel_own( p*p/4. - q );
    printf("Die Lösungen der quadratischen Gleichung %fx² + %fx +%f sind:\nx_1 = %f\nx_2 = %f\n", a, b, c, x1, x2);
}

int main(){
    double a = 5.5;
    double n = 2.0;
    double np = 2.0;
    double ap = 5.5;
    printf("Das Quadrat von %f ist %f.\n", n, square_to(&np));
    printf("Die Wurzel von %f ist %f.\n", a, root_to(&ap));
    printf("Und zuletzt schauen wir uns die Variablen an: %f und %f.\n", ap, np);
    square_equ(1., 10., 3.);
    printf("Die Wurzel von  0 ist %f\n", wurzel_own(0.));
    return 0;
}
