#include <stdio.h>
#include <math.h>

double root(double a){
    double a_i = a;
    double a_i1 = 0.5*(a_i + a/a_i);
    while(fabs(a_i1 - a_i) > 1e-10){
        a_i = a_i1;
        a_i1 = 0.5*(a_i + a/a_i);

    }
    return a_i1;
}
int main(){
    double m = 9.;
    printf("Die Wurzel von %f ist %f.\n", m, root(m));
    return 0;
}