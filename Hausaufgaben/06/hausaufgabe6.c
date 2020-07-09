// Dirk Knott und Lars Döpper
// make
// ./hausaufgabe6

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "filehelper.h"
#include "arrayhelpers.h"
#include "random.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double signal(double t, double t_max, double alpha){
    if( (t>=0) && (t<=t_max) ) {
        return sin(2.0*M_PI*alpha*t);
    }
    else {
        return 0.0;
    }
}

double echo(double t, double tl, double t_max, double alpha, double a, double b, double beta){
    double sign, rand, fak_sin;
    fak_sin = b*sin(2.0*M_PI*beta*t);
    rand = 2.0*a*(genrand_res53() - 0.5);
    sign = signal(t-tl, t_max, alpha);
    return sign+rand+fak_sin;
}

int main(int argc, char* argv[]){
    unsigned long seed = 123456;
    init_genrand(seed);

    double *t, *e;
    double amp = 0.25;
    int length = 10001;
    char* files[] = {"amp0_5.txt", "amp1_0.txt", "amp2_0.txt", "amp4_0.txt"};
    t = (double*)malloc(sizeof(double)*length);
    e = (double*)malloc(sizeof(double)*length);
    for(int j=0; j<4; j++){
        amp *=2.0;
        for(int i=0; i<length; i++){
            t[i] = i*0.01;
            e[i] = echo(t[i], 50.0, 5.0, 10.0, amp, amp, 1.0);
            e[i] *= e[i];
        }
        print_data_table_double(files[j], t, e, length);
    }
    free(t); free(e);
    return 0;
}