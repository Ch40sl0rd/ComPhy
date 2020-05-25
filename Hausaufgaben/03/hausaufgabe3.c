#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "arrayhelpers.h"
#include "numerik_own.h"
#include "filehelper.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct numerov_param{
    double *f_array;
    double *g_array;
    double *s_array;
    int steps;
    double h;
}numerov_param;

double s_func(double x){
    return 0.0;
}

double g_func(double x){
    return 0.0;
}

/*******************************************************************************************************************************
 *  This function creates all needed datapoints of the functions g(x) und s(x)
 *  to solve the differential equation:
 *  f(x)'' + g(x)*f(x) = s(x)
 *  
 *  \param start: starting point of the interval
 *  \param end: end point of the interval
 *  \param steps: number of data points
 *  \param g: pointer to an array to be filled with datapoints for g(x).
 *  \param g_func: function g(x)
 *  \param s: pointer to an array to be filled with datapoints for s(x).
 *  \param s_func: function s(x).
 * 
 *  \return step-width h.
 * *****************************************************************************************************************************/
double* numerov_init(double start, double end, int steps, double *g, double(*g_func)(double), double *s, double(*s_func)(double) ){
    double* x_array;
    x_array = create_array_double(start, end, steps);

    for(int i=0; i<steps; i++){
        g[i] = g_func(x_array[i]);
        s[i] = s_func(x_array[i]);
    }
    return x_array;
}

void numerov_up(numerov_param parameters, double f_0, double f_1){
    double fac_u_n, fac_u_nm1, fac_u_np1;
    double fac_s;
    double *f, *g, *s, h;
    int steps;

    f = parameters.f_array;
    g = parameters.g_array;
    s = parameters.s_array;
    h = parameters.h;
    steps = parameters.steps;

    f[0] = f_0;
    f[1] = f_1;
    for(int i=1; i<steps-1; i++){
        fac_u_np1 = 1.0 + h*h/12.0*g[i+1];
        fac_u_n = 2.0*(1.0 - 5*h*h/12.0*g[i]);
        fac_u_nm1 = 1.0 + h*h/12.0 * g[i-1];
        fac_s = h*h/12.0*(s[i+1]+ 10.0*s[i] + s[i-1]);
        f[i+1] = 1/fac_u_np1*(fac_s + fac_u_n*f[i] - fac_u_nm1*f[i-1]);
    }
}

void numerov_down(numerov_param parameters, double f_max, double f_maxm1){
    double fac_u_n, fac_u_nm1, fac_u_np1;
    double fac_s;
    double *f, *g, *s, h;
    int steps;

    f = parameters.f_array;
    g = parameters.g_array;
    s = parameters.s_array;
    h = parameters.h;
    steps = parameters.steps;

    f[steps-1] = f_max;
    f[steps-2] = f_maxm1;
    for(int i=steps-2; i<1; i--){
        fac_u_np1 = 1.0 + h*h/12.0*g[i+1];
        fac_u_n = 2.0*(1.0 - 5*h*h/12.0*g[i]);
        fac_u_nm1 = 1.0 + h*h/12.0 * g[i-1];
        fac_s = h*h/12.0*(s[i+1]+ 10.0*s[i] + s[i-1]);
        f[i-1] = 1/fac_u_nm1*(fac_s + fac_u_n*f[i] - fac_u_np1*f[i+1]);
    }
}

double bound_con_up(double free_param, numerov_param parameters, double f_0, double f_max){
    numerov_up(parameters, f_0, free_param);
    return parameters.f_array[parameters.steps] - f_max;
}

double bound_con_down(double free_param, numerov_param parameters, double f_0, double f_max){
    numerov_down(parameters, f_max, free_param);
    return parameters.f_array[0] - f_0;
}

/******************************************************************************************************************************************************************
 *  This function searches for a zero-crossing of a given numerov_function to determine the free parameter of the numerov-function.
 *  
 *  \param x0: first starting value for free parameter.
 *  \param x1: second starting value for free parameter.
 *  \param paramters: the parameters for the numerov function.
 *  \param func: numerov boundary condition function, either up or down
 *  \param f_0: value of f(start)
 *  \param f_max: value of f(end)
 *  \param steps: maxmimum number of steps for searching zero crossings
 
 * *****************************************************************************************************************************************************************/
double secant_numerov(double x0, double x1, numerov_param parameters, double(*func)(double, numerov_param , double, double), double f_0, double f_max, int steps){
    const double acc = 1e-12;
    double x_n;
    int step = 0;
    do{
        x_n = x1 - (x1-x0)*func(x1, parameters, f_0, f_max)/(func(x1, parameters, f_0, f_max) - func(x0, parameters, f_0, f_max));
        x0 = x1;
        x1 = x_n;
        steps++;
    }
    while(fabs(x0-x1)>acc && step<steps);
    return x1;
}
/*******************************************************************************************************
 *  This function exectues the numerov-methode to solve a differential euqation of the form
 *  f(x)'' + g(x)*f(x) = s(x) with the parameters f(start) = f_0; f(end) = f_max
 * 
 *  \param start: starting point of the interval
 *  \param end: endpoint of the interval
 *  \param steps: number of steps.
 *  \param g_func: pointer to function g(x)
 *  \param s_func: pointer to function s(x)
 *  \param f_0 value of f at start.
 *  \param f_max: value of f at end
 *  \param direction: either 1 for up or -1 for down

 *******************************************************************************************************/
double** numerov_complete(double start, double end , int steps, double (*g_func)(double), double(*s_func)(double), double f_0, double f_max, int direction){
    double *f_array, *g_array, *s_array, *x_array;
    double **data_table;
    double free_param, set_param;
    numerov_param paramters;
    void (*numerov_func)(numerov_param, double, double);
    double (*bound_con)(double, numerov_param, double, double);
    if(direction==1){
        numerov_func = numerov_up;
        bound_con = bound_con_up;
        set_param = f_0;
    }
    else if(direction == -1){
        numerov_func = numerov_down;
        bound_con = bound_con_down;
        set_param = f_max;
    }
    else{
        printf("[numerov_complete] No valid option. Function will abort. Either Choose 1 for forwards or -1 for backwards.\n");
        return NULL;
    }
    //create array for datapoints x_n and corresponding g(x) and s(x)
    g_array = (double*)malloc(sizeof(double)*steps);
    s_array = (double*)malloc(sizeof(double)*steps);
    x_array = numerov_init(start, end, steps, g_array, g_func, s_array, s_func);
    //allocate memory for values of f(x_n)
    f_array = (double*)malloc(sizeof(double)*steps);

    //create struct with all parameters
    paramters.f_array = f_array;
    paramters.g_array = g_array;
    paramters.s_array = s_array;
    paramters.steps = steps;
    paramters.h = (end-start)/(double)(steps-1);

    //solve differential equation for given boundary conditions
    free_param = secant_numerov(-1.0, 1.0, paramters, bound_con, f_0, f_max, 12);

    //use free parameter to calculate final values of the function
    numerov_func(paramters, set_param, free_param);
    data_table = create_2d_array(steps, 2);
    for(int i=0; i<steps; i++){
        data_table[i][0] = x_array[0];
        data_table[i][1] = f_array[0];
    }


    free(g_array); free(s_array); free(x_array); free(f_array);
    return data_table;
}

int main(int argc, char* argv[]){
    printf("Hallo Welt. test\n");
    return 0;
}