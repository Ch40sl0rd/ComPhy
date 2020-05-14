// make
// ./hausaufgabe1
// Lars Döpper
#include <stdio.h>
#include <gsl/gsl_integration.h>
#include "filehelper.h"
double derivate_sym_one(double x, double h, double(*func)(double)){
    return (func(x+h)-func(x-h))/(2*h);
}

double* derivate_sym_one_array(double* data, double a, double b, const int length){
    double* data_ret;
    double h;
    if(data == NULL){
        printf("Error while handling an empty pointer.\n");
        exit(404);
    }
    data_ret = (double*)malloc(sizeof(double)*length);
    if(!data_ret){
        printf("Error while allocating local memory.\n");
        exit(1);
    }
    h = (b-a)/(length-1);
    //use normal one-sided derivative for the first and last data point
    data_ret[0] = (data[1]-data[0])/h;
    data_ret[length-1] = (data[length-1]-data[length-2])/h;
    for(int i=1; i<length-1; i++){
        data_ret[i] = (data[i+1]-data[i-1])/(2.*h);
    }
    return data_ret;
}

double derivate_sym_two(double x, double h, double(*func)(double)){
    return (func(x+h) + func(x-h) -2*func(x))/(h*h);
}

double* derivate_sym_two_array(double* data, double a, double b, const int length){
    double* data_ret;
    double h;
    if(data == NULL){
        printf("Error while handling an empty pointer.\n");
        exit(404);
    }
    data_ret = (double*)malloc(sizeof(double)*length);
    if(!data_ret){
        printf("Error while allocating local memory.\n");
        exit(1);
    }
    h = (b-a)/(length-1);
    //use normal one-sided derivative for the first and last data point
    data_ret[0] = (data[2]-2.*data[1]+data[0])/(h*h);
    data_ret[length-1] = (data[length-1]- 2.*data[length-2] + data[length-3])/(h*h);
    for(int i=1; i<length-1; i++){
        data_ret[i] = (data[i+1] + data[i-1] -2.*data[i])/(h*h);
    }
    return data_ret;
}

double derivate_sym_three(double x, double h, double(*func)(double)){
    return ( (func(x+2*h) - func(x-2*h)) -2*(func(x+h) - func(x-h)) )/(2*h*h*h);
}

double zero_crossing(double (*func)(double, void*), const double x0, const double x1, double acc, void *p){
    int max_steps = 100;
    double x_i, x_ip1, x_im1;
    int steps = 0;
    x_i = x1;
    x_im1 = x0;
    do{
        x_ip1 = x_i - func(x_i, p)* (x_i - x_im1)/(func(x_i, p) - func(x_im1, p));
        x_im1 = x_i;
        x_i = x_ip1;
        steps++;
        //printf("Momentane geschätzte Position der Nullstelle: %15.6e\n", x_ip1);
    }
    while( (fabs(x_ip1 - x_im1) > acc) && steps < max_steps);
    if(steps ==max_steps){
        printf("Maximum number of steps has been reached. No zero crossing detected.\n");
        return 0.;
    }
    //printf("This function needed %d steps of iteration.\n", steps);
    return x_ip1;
}

double integrate_trapez(double (*f)(double, void*), double a, double b, int n, void *p){
    double h;
    h = (b-a)/(double)(n-1);
    double sum;
    sum = h/2*(f(a, p)+f(b, p));
    for(int i=1; i<n-1; i++){
        sum += h*f(a+i*h, p);
    }
    return sum;
}

double integrate_trapez_adap( double(*f)(double, void *), double a, double b, double acc, void *p){
    double sum_i, sum_ip1;  //zwei aufeinanderfolgende integralsummen
    double h;   //Schrittweite
    int N; //Anzahl der Schritte
    
    N = 10; //Startschrittzahl
    h = fabs(b-a)/(double)(N-1); //startschrittweite
    do{
        sum_i = sum_ip1;
        sum_ip1 = h/2*(f(a, p) + f(b, p));
        for(int i=1; i<N; i++){
            sum_ip1 += h*f(a+i*h, p);
        }

        N *= 2;
        h = (b-a)/(double)(N-1);
        //printf("Der Momentane Wert des Integrals ist %15.6e \n", sum_ip1);
    }
    while(fabs(sum_i - sum_ip1) > acc);
    return sum_ip1;
}

int index_rom(int i, int k, int m_max){
    return (i-k + k*(m_max+1) - k*(k-1)/2);
}

double integrate_romberg(double a, double b, int n_0, double (*func)(double, void*), double acc, int m_max, void* p){
    int k, m, n, m_now; //Indezies
    double *h; //list of the step width in each iteration.
    double *T_int; //list of the values of the trapez integration.
    double *T_tilde; //list for t_tilde function after the Neville-scheme.
    double res; //value of the integrale.

    //allocating all needed memory.
    h = (double*)malloc( (m_max+1)*sizeof(double));
    T_int = (double*)malloc((m_max+1)*sizeof(double));
    T_tilde = (double*)malloc(sizeof(double)* (m_max+1)*(m_max+2)/2);

    if(!h || !T_int || !T_tilde){
        printf("Error while allocating memory.\n");
        exit(1);
    }

    //start with first step of iteration:
    n = n_0;
    h[0] = (b-a)/(double)(n-1);
    T_int[0] = integrate_trapez(func, a, b, n, p);
    T_tilde[index_rom(0,0, m_max)] = T_int[0];
    //calculate all other T_tilde accoreding to the Neville-scheme
    for(m=1; m<=m_max; ++m){
        //calculate number of steps and step width.
        n *= 2;
        h[m] = (b-a)/(double)(n-1);
        //calculate new value of the trapez integrale.
        T_int[m] = integrate_trapez(func, a, b, n, p);
        T_tilde[index_rom(m,0, m_max)] = T_int[m];

        //generate following T_tilde after the Neville-scheme
        for(k=1; k<=m; k++){
            T_tilde[index_rom(m,k, m_max)] = -(h[m-k]*h[m-k])/(h[m]*h[m]-(h[m-k]*h[m-k])) * T_tilde[index_rom(m, k-1, m_max)]
                                         + (h[m]*h[m])/(h[m]*h[m]-(h[m-k]*h[m-k])) * T_tilde[index_rom(m-1, k-1, m_max)];
        }
        //printf("%5d\t%15.6e\t%15.6e\n", m, T_tilde[index_rom(m,m, m_max)], T_int[m]);
        m_now = m;
        //end precedure if accuracy is reached.
        if(fabs(T_tilde[index_rom(m,m,m_max)] - T_tilde[index_rom(m-1, m-1, m_max)]) <= acc){
            break;
        }
        
        
    }
    if(m>m_max){
        printf("No convergence reached in %d steps\n", m);
        free(T_tilde); free(T_int); free(h);
        return NAN;
    }
    res = T_tilde[index_rom(m_now,m_now,m_max)];
    //printf("Result reached after %d steps\n", m);
    free(T_tilde); free(T_int); free(h);
    return res;
}

double  gamma_integrand(double t, void *p){
    double z = ((double*)p)[0];
    return exp(-t)*pow(t, z-1);
}

double gamma_func(double z){
    if( z>2){
        return  (z-1)*gamma_func(z-1);
    }
    double a, b, acc;
    a = 0.;
    b = 10*z;
    acc = 1e-8;
    return (integrate_trapez_adap(&gamma_integrand, a, b, acc, (void*)(&z)));
}

void gaus_legendre(double a, double b, double *x_vals, double *w_vals, int n){
    gsl_integration_glfixed_table *xw_table;
    size_t i;

    xw_table = gsl_integration_glfixed_table_alloc(n);
    if(!xw_table){
        printf("Error while handling Gauss-Legendre integration.\n");
        exit(4);
    }

    for(i =0; i<n; i++){
        gsl_integration_glfixed_point(a, b, i, &x_vals[i], &w_vals[i], xw_table);
    }
    gsl_integration_glfixed_table_free(xw_table);
}

double gl_integrate(double a, double b, double *x, double *w, int n, double(*f)(double, void *), void *p){
    double sum;
    gaus_legendre(a,b, x,w, n);
    sum = 0.;
    for(int i=0; i<n; i++){
        sum += f(x[i], p)*w[i];
    }
    return sum;
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

/********************************************************************************************************************
 *  This function implements on step of the runge-kutta-methode of the third order.
 * 
 *  \param neq: order of differential equation (deq)
 *  \param h: width of one step.
 *  \param t: current time.
 *  \param y: vector of values of y_n at point t 
 *  \param f: vector of values for f_n at time t.
 *  \param p: list of parameters for deq-function
 *  \param dgl_func: function that represents the differential equation. 
 * 
 * ********************************************************************************************************************/
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
        val_array[i][0] = t;
        for(int j=1; j<=neq; j++){
            val_array[i][j] = y[j-1];
        }
    }
    print_table2file(NULL, val_array, neq+1, step_number);
    free(f); free(y); free(values); free(val_array);
}