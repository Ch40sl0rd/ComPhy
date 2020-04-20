#include <stdio.h>
#include <math.h>

int main(){
    double k = 2.;
    double sum_k =  1.;
    while(fabs(1/(k*k)) > 1e-10){
        sum_k += 1/(k*k);
        k++;
    }
    printf("Der Wert der Reihe  ist %f.\n", sum_k);
    printf("Der theoretische Wert ist %f.\n", (3.1415*3.1415)/6);
    return 0;
}