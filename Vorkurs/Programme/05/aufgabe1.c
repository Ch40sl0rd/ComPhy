#include <stdio.h>
#include <math.h>
#include "arrayhelpers.h"
typedef int myint;

unsigned long int random(int I_j){
    const int m = 2147483647;
    const int a = 16807;
    const int q = 127773;
    const int r = 2836;
    int I_jp1 = a*(I_j%q) - r*(I_j/q);
    if(I_jp1 <= 0){
        I_jp1 += m;
    }
    return  I_jp1;
}

double rand_frac(int d){
    return (double)d/(double)2147483647;
}

int main(){
    int  d =  156874;
    const int n = 1000;
    for(unsigned int i = 0; i<n; i++){
        d = random(d);
        printf("%d %f\n", d, rand_frac(d));
        
    }
}