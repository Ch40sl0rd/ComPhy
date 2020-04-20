#include <stdio.h>
#include <math.h>

long long fakultaet(int n){
    long  long ergebnis = 1;
    while(n>1){
        ergebnis*=n;
        n--;
    }
    return ergebnis;
}

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

double power_own(double x, int n){
    if(n==0){
        return 1;
    }
    if(n%2==0){
        return power_own(x*x, n/2);
    }
    else{
        return x*power_own(x*x, (n-1)/2);
    }
}

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

double zeta(double s){
    double ret=0.;
    double part_sum = 1.;
    for(int i=1; part_sum>1e-10;i++){
        part_sum = pow( (double)i, -s);
        ret += part_sum;
    }
    return ret;
}