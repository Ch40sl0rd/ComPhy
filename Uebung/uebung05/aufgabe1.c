#include <stdio.h>
#include <stdlib.h>

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

void euler_step(int neq, double h, double t, double *y, double *f, void *p,
    void (*dgl_func)(int, double, double*, double*, void*))
    {
        //führe Euler-schritt durch
        dgl_func(neq, t, y, f, p);

        //calculate new values for y_n+1
        for(int i=0; i<neq; i++){
            y[i] += h*f[i];
        }
    }

void runge_kutta_2_step(int neq, double h, double t, double *y, double *f, void *p,
    void (*dgl_func)(int, double, double*, double*, void*))
    {
    //calculate the new value of y_np1 following the rk2 metode
    double *k;
    k = (double*)malloc(sizeof(double)*neq);
    //save values for
    dgl_func(neq, t, y, f, p);
    for(int i=0;i<neq;i++){
        k[i] = y[i] + h/2*f[i];
    }
    dgl_func(neq, t+0.5*h, k, f, p);
    for(int i=0; i<neq; i++){
        y[i] = y[i] + h*f[i];
    }
    free(k);
}

void runge_kutta_3_step(int neq, double h, double t, double *y, double *f, void *p,
    void (*dgl_func)(int, double, double*, double*, void*))
{
    double *k1, *k2, *k3, *yh;
    k1 = (double*)malloc(sizeof(double)*neq);   
    k2 = (double*)malloc(sizeof(double)*neq);
    k3 = (double*)malloc(sizeof(double)*neq);
    yh = (double*)malloc(sizeof(double)*neq);   //help field for storing values of y_n
    //save value for f_n in f by calculating vlaue of dgl_function
    dgl_func(neq, t, y, f, p);
    for(int i=0; i<neq; i++){
        k1[i] = h*f[i];
        yh[i] = y[i] + 0.5*k1[i];
    }
    dgl_func(neq, t+0.5*h, yh, f, p);
    for(int i=0; i<neq; i++){
        k2[i] = h*f[i];
        yh[i] = y[i]-k1[i] + 2*k2[i];
    }
    dgl_func(neq, t+h, yh, f, p);
    for(int i=0; i<neq; i++){
        k3[i] = h*f[i];
    }
    for(int i=0; i<neq; i++){
        y[i] += (k1[i] + 4*k2[i] + k3[i])/6.0;
    }
    free(k1); free(k2); free(k3); free(yh);
}

void runge_kutta_4_step(int neq, double h, double t, double *y, double *f, void *p,
    void (*dgl_func)(int, double, double*, double*, void*))
{
    double *k1, *k2, *k3, *k4, *yh;
    k1 = (double*)malloc(sizeof(double)*neq);
    k2 = (double*)malloc(sizeof(double)*neq);
    k3 = (double*)malloc(sizeof(double)*neq);
    k4 = (double*)malloc(sizeof(double)*neq);
    yh = (double*)malloc(sizeof(double)*neq);

    dgl_func(neq, t, y, f, p);
    for(int i=0; i<neq; i++){
        k1[i] = h*f[i];
        yh[i] = y[i]+k1[i]/2;
    }
    dgl_func(neq, t+h/2.0, yh, f, p);

    for(int i=0; i<neq; i++){
        k2[i] = h*f[i];
        yh[i] = y[i] + k2[i]/2.0;
    }
    dgl_func(neq, t+h/2.0, yh, f, p);
    for(int i=0; i<neq; i++){
        k3[i] = h*f[i];
        yh[i] = y[i] + k3[i];
    }
    dgl_func(neq, t+h, yh, f, p);
    for(int i=0; i<neq; i++){
        k4[i] = h*f[i];
    }

    for(int i=0; i<neq; i++){
        y[i] += (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.0;
    }

    free(k1);free(k2);free(k3);free(k4);free(yh);
}

/*******************************************************************************************************************
 *  This function solves a given differantial euqation in vector noation by apllying either
 *  the euler-cauchy-methode or a runge-kutta-methode of different orders to solve the problem.
 *  Starting values have to given as well.
 *  This function gives out all vales to the console or file(not implemented yet).
 * 
 *  \param neq: Order of deq
 *  \param start \param stop start- and endpoint of the unknown function.
 *  \param step_size: step size between two points.
 *  \param number_methode: Number of the methode to be apllied. 1 for euler.cauchy or 2-4 for rk-methode
 *  \param p: List of all parameters for the deq-function.
 *  \param start_vals: List of starting values for the unknown function.
 *  \param dgl_func: vector-valued function of the differantiel equation 
 * 
 ***********************************************************************************************************************/
void solve_dgl(int neq, double start, double stop, double step_size, int number_methode, void *p, double *start_vals,
                void (*dgl_func)(int, double, double*,double*, void*))
{
    int step_number;
    double *y, *f,  t, *values, **val_array;
    void (*dgl_solve_step)(int, double, double, double*, double*,void*, void (*func)(int, double, double*, double*, void*));

    switch(number_methode){
        case 1: printf("Euler-Verfahren"); dgl_solve_step = &euler_step; break;
        case 2: printf("Runge-Kutta-Verfahren der Stufe 2.\n");dgl_solve_step=&runge_kutta_2_step; break;
        case 3: printf("Runge-Kutta-Verfahren der Stufe 3.\n");dgl_solve_step=&runge_kutta_3_step; break;
        case 4: printf("Runge-Kutta-Verfahren der Stufe 4.\n");dgl_solve_step=&runge_kutta_4_step; break;
        default: printf("Please chose an integrated methode to solve the deq.\n"); abort(); break;
    } 

    step_number = (stop-start)/step_size +1;

    //allocate 2d array to save values for each step
    values = (double*)malloc(sizeof(double*)*(neq+1)*step_number);
    val_array = (double**)malloc(sizeof(double*)*step_number);
    for(int i=0; i<step_number; i++){
        val_array[i] = values+i*(neq+1);
    }
    //allocate memory for vectors for y and f
    y = (double*)malloc(sizeof(double)*neq);
    f = (double*)malloc(sizeof(double)*neq);

    //start-values
    for(int i=0; i<neq; i++){
        y[i] = start_vals[i];
    }
    t = start;
    //print out starting values
    printf("%15.6e\t%15.6e\t%15.6e\n", t, y[0], y[1]);

    for(int i=1; i<step_number; i++){
        dgl_solve_step(neq, step_size, t, y, f,p, dgl_func);
        t += step_size;
        printf("%15.6e\t%15.6e\t%15.6e\n", t, y[0], y[1]);
        val_array[i][0] = t;
        for(int j=1; j<=neq; j++){
            val_array[i][j] = y[j-1];
        }
    }
    free(f); free(y); free(values); free(val_array);
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