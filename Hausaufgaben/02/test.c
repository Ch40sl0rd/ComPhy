#include <stdio.h>
#include <stdlib.h>
#include "filehelper.h"
#include "numerik_own.h"
#include "h2.h"

/**
 * This function calculates the value of the function of the 
 * differential equation to be used in solving differential equation
 * This is a function for an harmonic oszilator.

 * \param neq:  degree of the differential equation.
 * \param t:    current time step
 * \param y:    list of values of vector y at the current time step
 * \param f:    list for calculating the value value of the function for
 *              the differential equation.
**/
void dgl_function(int neq, double t, double *y, double *f, void *p){
    if(neq>2){
        printf("Falscher Grad der DGL.\n");
        abort();
    }
    
    f[0] = y[1];
    f[1] = -y[0];
}

int main(int argc, char** argv){
    double start, stop, h, *start_vals;
    int neq;

    start = 0.;
    stop = 3.1415;
    h = 1e-3;
    neq = 2;
    start_vals = (double*)malloc(sizeof(double)*neq);

    if(argc<2){
        start_vals[0] = 0.0;
        start_vals[1] = 1.0;
    }
    else{
        start_vals[0] = atof(argv[1]);
        start_vals[1] = atof(argv[2]);
    }

    solve_dgl(neq, start, stop, h, 4,NULL,start_vals,  &dgl_function);
    free(start_vals);
}