#include <stdio.h>
#include <math.h>

int sgn_own(double x){
    if(x<0){
        return -1;
    }
    else if(x==0){
        return 0;
    }
    else{
        return 1;
    }

}

double betrag_own(double x){
    if(x==0){
        return 0.;
    }
    else if(x>0){
        return x;
    }
    else{
        return -1*x;
    }
}

double cos_own(double x){
    double part_sum = 1.;
    double ergebnis = 1.;
    int sign = 1;
    for(int i = 1; i <  100; i++){
        sign*=-1;
        part_sum *= x*x/((2.*i-1)*2.*i);
        ergebnis+= sign*part_sum;
    }
    return ergebnis;
}

double wurzel_own(double x){
    double a_i = x;
    double a_i1 = 0.5*(a_i + x/a_i);
    while(fabs(a_i1 - a_i) > 1e-10){
        a_i = a_i1;
        a_i1 = 0.5*(a_i + x/a_i);

    }
    return a_i1;
}

int main(){
    double cos_test = 3.1415;
    double sgn_test = -0.5;
    double betrag_test = -0.1;
    double wurzel_test = 1586542.;
    printf("Test von Cosinus mit %f = %f .\n", cos_test, cos_own(cos_test));
    printf("Test von Sign mit %f = %d .\n", sgn_test, sgn_own(sgn_test));
    printf("Test von Betrag mit %f = %f .\n", betrag_test, betrag_own(betrag_test));
    printf("Test von Wurzel mit %f = %f .\n", wurzel_test, wurzel_own(wurzel_test));
}