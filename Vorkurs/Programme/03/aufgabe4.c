#include <stdio.h>
#include <math.h>

double exp_own(double x){
    double part_sum = 1.;
    double ret = 1.;
    for(int i=1; part_sum>1e-11;i++){
        part_sum *= x/i;
        ret += part_sum;
    }
    return ret;
}

double log_own(double x){
    if(x<=0.){
        printf("Dieser Wert ist nicht zulässig\n");
        return -1;
    }
    double ai = (x-1)/(x+1);
    double ret = ai;
    //ai+1  = ai*(x-1)/(x+1)**2*i/i+2
    for(int i = 1; ai>1e-11; i++){
        ai *= pow((x-1)/(x+1), 2.)*i/(i+2);
        ret += ai;
    }
    return 2*ret;
}

double  power(double x, double y){
    double z = log_own(x);
    return exp_own(y*z);
}

int main(){
    double x  =  1.;
    double y = 1.;
    printf("Test  der power-func  %f\n", power(x,y));
}