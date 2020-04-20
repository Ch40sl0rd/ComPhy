#include <stdio.h>
#include <math.h>

double zeta(double s){
    double ret=0.;
    double part_sum = 1.;
    for(int i=1; part_sum>1e-10;i++){
        part_sum = pow( (double)i, -s);
        ret += part_sum;
    }
    return ret;
}

int main(){
    double x = 1.5;
    printf("test der Ceta-Funktion %f\n", zeta(x));
}