//Lars Döpper, Dirk Knott
//make
// ./hausaufgabe3
//Das Programm benötigt allerdings sehr viel Arbeitsspeicher
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "arrayhelpers.h"
#include "numerik_own.h"
#include "filehelper.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double lambda; //paramters for calculating eigenvalues
double epsilon; //paramters for calculating electric potential

//define struct for paramters for the numerov-methode to shorten
//function calls
typedef struct numerov_param{
    double *f_array;
    double *g_array;
    double *s_array;
    int steps;
    double h;
}numerov_param;

double s_const_0(double x){
    return 0.0;
}

double eigenvalue_function(double x){
    return -lambda;
}

double pot_period_electric(double x){
    return -2.0*(60*pow(cos(M_PI*x), 16) - lambda + epsilon*x);
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
    x_array = (double*)malloc(sizeof(double)*steps);
    double h = (end-start)/(steps-1);

    for(int i=0; i<steps; i++){
        x_array[i] = start + i*h;
        g[i] = g_func(x_array[i]);
        s[i] = s_func(x_array[i]);
    }
    return x_array;
}

/*************************************************************************************
 *  This function uses numerov's methode to calculate the value of the function f
 *  which solves the differential equation:
 *  f(x)'' + g(x)f(x) = s(x)
 *  This methode starts at the first two values and goes forwards.
 * 
 *  \param paramters: paramters to be used by other numerov functions
 *  \param f_max: value of f at the first point
 *  \param f_maxm1: value of f at the second point.
 *************************************************************************************/
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
        f[i+1] = (fac_s + fac_u_n*f[i] - fac_u_nm1*f[i-1])/fac_u_np1;
    }
}

/*************************************************************************************
 *  This function uses numerov's methode to calculate the value of the function f
 *  which solves the differential equation:
 *  f(x)'' + g(x)f(x) = s(x)
 *  This methode starts at the last two values and goes backwards.
 * 
 *  \param paramters: paramters to be used by other numerov functions
 *  \param f_max: value of f at the last point
 *  \param f_maxm1: value of f at the penultimate point.
 *************************************************************************************/
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
        f[i-1] = (fac_s + fac_u_n*f[i] - fac_u_np1*f[i+1])/fac_u_nm1;
    }
}

/*******************************************************************************************
 *  This function calculates the difference between the current value at the last point
 *  and the set value at this point. It is mainly used by numerov_complete to find the ideal
 *  free paramter.
 * 
 *  \param free_param: the free parameter at second position to define the ideal function
 *  \param paramters: parameters for the other numerov_methods
 *  \param f_0: boundary condition f(start)
 *  \param f_max: boundary condition f(end)
 * 
 *  \return difference betwwen f(start) and f_max
 * *****************************************************************************************/
double bound_con_up(double free_param, numerov_param parameters, double f_0, double f_max){
    numerov_up(parameters, f_0, free_param);
    //printf("Die Differenz beträgt momentan %15.6e\n", parameters.f_array[parameters.steps-1] - f_max);
    return parameters.f_array[parameters.steps-1] - f_max;
}

/*******************************************************************************************
 *  This function calculates the difference between the current value at the first point
 *  and the set value at this point. It is mainly used by numerov_complete to find the ideal
 *  free paramter.
 * 
 *  \param free_param: the free parameter at second position to define the ideal function
 *  \param paramters: parameters for the other numerov_methods
 *  \param f_0: boundary condition f(start)
 *  \param f_max: boundary condition f(end)
 * 
 *  \return difference betwwen f(start) and f_max
 * *****************************************************************************************/
double bound_con_down(double free_param, numerov_param parameters, double f_0, double f_max){
    numerov_down(parameters, f_max, free_param);
    return (parameters.f_array[0] - f_0);
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
    double x_n, x_np1, temp;
    double f_xn, f_xnp1;
    int step = 0;
    x_n = x0;
    x_np1 = x1;
    do{
        //temp = x_np1 - (x_np1-x_n)/(func(x_np1, parameters, f_0, f_max) - func(x_n, parameters, f_0, f_max))*func(x_np1, parameters, f_0, f_max);
        f_xn = func(x_n, parameters, f_0, f_max);
        f_xnp1 = func(x_np1, parameters, f_0, f_max);
        //printf("Die Werte der Funktion betragen jeweils: %15.6e; %15.6e\n", f_xnp1, f_xn);
        temp = x_np1 - f_xnp1*(x_np1 - x_n)/(f_xnp1 - f_xn);
        //temp = x_np1-func(x_np1, parameters, f_0, f_max)*(x_np1-x_n)/(func(x_np1, parameters, f_0, f_max) - func(x_n, parameters, f_0, f_max));
        //printf("Der neue Schätzwert lautet: %15.12e\n", temp);
        x_n = x_np1;
        x_np1 = temp;
        step++;
    }
    while(fabs(x_np1 - x_n)>acc && step<steps);
    //printf("Konvergenz erreicht nach %d Schritten.\n", step);
    return x_np1;
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
    numerov_param parameters;
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
    /**print_data2file(NULL, x_array, steps);
    print_data2file(NULL, g_array, steps);
    print_data2file(NULL, s_array, steps);
    //allocate memory for values of f(x_n)**/
    f_array = (double*)malloc(sizeof(double)*steps);

    //create struct with all parameters
    parameters.f_array = f_array;
    parameters.g_array = g_array;
    parameters.s_array = s_array;
    parameters.steps = steps;
    parameters.h = (end-start)/(double)(steps-1);

    //print_data2file(NULL, parameters.g_array, steps);
    //print_data2file(NULL, parameters.s_array, steps);

    //solve differential equation for given boundary conditions
    free_param = secant_numerov(-1.0, 1.5, parameters, bound_con, f_0, f_max, 20);
    //printf("[numerov_complete] Der freie Parameter lautet : %15.6e\n", free_param);

    //use free parameter to calculate final values of the function
    numerov_func(parameters, set_param, free_param);
    //print_data2file(NULL, x_array, steps);
    data_table = create_2d_array(steps, 2);
    for(int i=0; i<steps; i++){
        data_table[i][0] = x_array[i];
        data_table[i][1] = f_array[i];
    }


    free(g_array); free(s_array); free(x_array); free(f_array);
    return data_table;
}

/**********************************************************************************************************************************
 *  This function searches for eigenvalues of the following differential equation:
 *  f''(x) + g(x)f(x) - lambdaf(x) = s(x) with set boundary conditions.
 *  it uses the numerov-methode in conjunktion with the secant-methode
 *  to find the coreesponding values for the function x and evaluates the eigenvalues
 *  from the basis of their maximum amplitude.
 *  
 *  \param start: starting point of the interval of eigenvalues
 *  \param end: end point of interval fro eigenvalues
 *  \param step_width: width between two possible eigenvalues
 *  \param g_function: function g(x) for the differential equation.
 *  \param s_function: function s(x) for differential equation.
 *  \param start_f: starting point for x
 *  \param end_f: end point for x
 *  \param f_steps: number of steps for x between start_f and end_f
 *  \param f_0: value for f(start)
 *  \param f_max: value for f(end)
 *  \param direction: either 1 or -1 for direction of numerov methode
 *  \param max_num_eigenvalues: number of eigenvalues to look for.
 * 
 *  \return list of all found eigenvalues.
 ************************************************************************************************************************************/
double* search_eigenvalues(double start, double end, double step_width, 
    double (*g_function)(double), double(*s_function)(double),
    double start_f, double end_f, int f_steps, double f_0, double f_end, int direction,
    int max_num_eigenvals){

    double *list_eigenvals, max, **function_values;
    int num_eigenval, i, num_steps;

    //set and calucalte constants for methode:
    double threshhold = 0.5e3;    //threshhold for maxmimum amplitude, al non-eigenvalues have around O(10) max amplitude
    num_steps = abs((int)((end-start)/step_width +1));
    //printf("die Anzahl der Schritte beträgt %d\n", num_steps);
    num_eigenval = 0;

    //allocate memory for list of eigenvalues
    list_eigenvals = (double*)malloc(sizeof(double)*max_num_eigenvals);

    for(i=0; i<num_steps; i++){
        lambda = start + i*step_width;
        //printf("Scan for lambda = %15.6e\n", lambda);
        function_values = numerov_complete(start_f, end_f, f_steps, g_function, s_function, f_0, f_end, 1);
        //print_table2file("test.txt", function_values, f_steps, 2);
        //find maxmimum amplitude of the function
        max = 0.0;
        for(int j=0; j<f_steps; j++){
            if(max < fabs(function_values[j][1])){
                max = fabs(function_values[j][1]);
            }
        }
        //printf("Das Maximum beträgt %15.6e.\n", max);
        if(max > threshhold){
            //printf("[search_eigenvals] The maxmimum is %15.6e\n", max);
            printf("[search_eigenvals] The %dth eigenvalue is %15.6e\n",num_eigenval+1, lambda);
            list_eigenvals[num_eigenval]=lambda;
            num_eigenval++;
            //increase i by a significant ammount so that the next chechek value is not the same eigenvalue
            i += 10;
        }
        free2d(function_values);
        //break the loop if the desired number of eigenvalues has been found.
        if(num_eigenval == max_num_eigenvals){
            i=num_steps;   
        }
    }
    printf("[search_eigenvals] A total of %d eigenvalues have been found.\n", num_eigenval);
    if(num_eigenval<max_num_eigenvals){
        for(int j=num_eigenval; j<max_num_eigenvals; j++){
            list_eigenvals[j] = 0.0;
        }
    }
    return list_eigenvals;
}

void plot_function_gap(double lambda_m, char* const filename){
    double **function_values;
    double sum;
    double h = 60.0/(1000.0-1.0);
    epsilon = 0.0;
    lambda = lambda_m;
    function_values = numerov_complete(0.0, 8.0, 1000, pot_period_electric, s_const_0, 0.0, 0.0, 1);
    //calculate norm-factor
    sum = 0.0;
    for(int i=0; i<1000-1; i++){
        sum += h/2*function_values[i][1]*function_values[i][1];
    }
    sum = sqrt(sum);
    for(int i=0; i<1000; i++){
        function_values[i][1] = function_values[i][1]/sum;
    }
    print_table2file(filename, function_values, 1000, 2);
    free2d(function_values);
}

double* search_eigenvals_pot(double start, double end, double width, double L, int* num_eigenvals){
    double* list_eigenvals;
    int cur_eigenval, max_num_eigenvals;
    int steps_eigenvals;
    double *g_array, *s_array, *x_array, *f_array;
    numerov_param parameters;
    double end_j, end_jm1;

    max_num_eigenvals = 300;
    cur_eigenval = 0;
    list_eigenvals = (double*)malloc(sizeof(double)*max_num_eigenvals);
    steps_eigenvals = abs((int)((end-start)/width)) +1;
    lambda = 0.0;

    f_array = (double*)malloc(sizeof(double)*1000);
    g_array = (double*)malloc(sizeof(double)*1000);
    s_array = (double*)malloc(sizeof(double)*1000);

    parameters.f_array = f_array;
    parameters.s_array = s_array;
    parameters.steps = 1000;
    parameters.g_array = g_array;
    parameters.h = fabs(L/(1000-1));
    
    end_j = 0;
    end_jm1 = 0;
    for(int i=0; i<steps_eigenvals; i++){
        lambda = start +i*width;
        x_array = numerov_init(0.0, L, 1000, g_array, pot_period_electric, s_array, s_const_0);
        numerov_up(parameters, 0.0, 0.000001);

        end_jm1 = end_j;
        end_j = parameters.f_array[999];
        if(end_j==0 || (end_j>0 && end_jm1<0) ||(end_j<0 && end_jm1>0) ){
            list_eigenvals[cur_eigenval] = lambda;
            cur_eigenval++;
            if(cur_eigenval==max_num_eigenvals){
                printf("Maximium number of eigenvalues have been reached.\n");
            }
        }
        free(x_array);
    }
    printf("[search_eigenvals_pot] A total of %d eigenvalues have been found.\n", cur_eigenval);
    *num_eigenvals=cur_eigenval;
    free(f_array); free(g_array); free(s_array);
    return list_eigenvals;
}



int main(int argc, char* argv[]){
    double *list_eigenvals_1, *list_eigenvals_2, *list_eigenvals_3, *list_eigenvals_4;
    int num_eigenvals;

    //##########################        Aufgabe 7.2     #######################
    printf("[main] Excersise 7.2 will be executed next, we search for the 10 lowest eigenvalues sorted by absolute value.\n");
    list_eigenvals_1 = search_eigenvalues(0.0, -0.3, -0.00001, eigenvalue_function, s_const_0, 0.0, 60.0, 1000, 1.0, 0.0, 1, 10);
    print_data2file("eigenvalues_1.txt", list_eigenvals_1, 10);
    free(list_eigenvals_1);

    //##########################        Aufgabe 8.1     #######################
    epsilon=0.0;
    printf("[main] Now calculating eigenvalues for periodic potential and L=8.\n");
    list_eigenvals_2 = search_eigenvals_pot(0.0, 60.0, 0.001, 8.0, &num_eigenvals);
    print_data2file("eigenvalues_l8.txt", list_eigenvals_2, num_eigenvals);
    free(list_eigenvals_2);
    plot_function_gap(5.960000e+00, "Gap1_lower.txt");
    plot_function_gap(1.926700e+01, "Gap1_higher.txt");
    plot_function_gap(2.339700e+01, "Gap2_lower.txt");
    plot_function_gap(3.974000e+01, "Gap2_higher.txt");

    //##########################        Aufgabe 8.3     #######################
    printf("[main] Now calculating eigenvalues for periodic potential and L variable.\n");
    list_eigenvals_3 = search_eigenvals_pot(0.0, 60.0, 0.001, 16.0, &num_eigenvals);
    print_data2file("eigenvalues_l16.txt", list_eigenvals_3, num_eigenvals);
    free(list_eigenvals_3);
    list_eigenvals_3 = search_eigenvals_pot(0.0, 60.0, 0.001, 32.0, &num_eigenvals);
    print_data2file("eigenvalues_l32.txt", list_eigenvals_3, num_eigenvals);
    free(list_eigenvals_3);
    list_eigenvals_3 = search_eigenvals_pot(0.0, 60.0, 0.001, 48.0, &num_eigenvals);
    print_data2file("eigenvalues_l48.txt", list_eigenvals_3, num_eigenvals);
    free(list_eigenvals_3);

    //##########################        Aufgabe 8.3     #######################
    printf("[main] Now calculating eigenvalues for electric potential\n");
    double eps_width = 0.2;
    char* list_txt[] = { "0_2.txt", "0_4.txt", "0_6.txt", "0_8.txt", "1_0.txt", "1_2.txt", "1_4.txt", "1_6.txt", "1_8.txt", "2_0.txt", };
    for(int i=1; i<=10; i++){
        epsilon = i*eps_width;
        list_eigenvals_4 = search_eigenvals_pot(0.0, 60.0, 0.001, 8.0, &num_eigenvals);
        print_data2file(list_txt[i-1], list_eigenvals_4, num_eigenvals);
        free(list_eigenvals_4);
    }
    return 0;
}