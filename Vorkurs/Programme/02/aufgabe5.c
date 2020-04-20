#include <stdio.h>
#include  <math.h>
#include <stdlib.h>

double cosin(double winkel, int n){
    double  cos_x  = 1.;
    double  x_sqr  = winkel *  winkel;
    int sign = 1;
    double part_sum = 1.;
    for(int i =1; i<n; i++){
        sign *= -1;
        part_sum  *= x_sqr /( (2*i-1)*2*i);
        cos_x += sign*part_sum;
    }
    return cos_x;
}
int main(int argc, char *argv[]){
    if(argc <  2){
        printf("Biite gebe einen Winkel  und einen Genauigkeitsgrad  an.\n");
        return -1;
    }
    double x = atof(argv[1]);
    int n = atoi(argv[2]);
    printf("Der Cosinus von %f aud math ist %f.\n",  x,  cos(x));
    printf("Der Cosinus von %f ist ungefähr %f.\n",  x, cosin(x,n));
    return 0;
}