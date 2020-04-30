#include <stdio.h>
#include <gsl/gsl_integration.h>
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
        printf("Momentane gewschätzte Position der Nullstelle: %15.6e\n", x_ip1);
    }
    while( (fabs(x_ip1 - x_im1) > acc) && steps < max_steps);
    if(steps ==max_steps){
        printf("Maximum number of steps has been reached. No zero crossing detected.\n");
        return 0.;
    }
    printf("This function needed %d steps of iteration.\n", steps);
    return x_ip1;
}

double integrate_trapez(double (*f)(double, void*), double a, double b, double h, void *p){
    double sum;
    sum = h/2*(f(a, p)+f(b, p));
    for(int i=1; a+i*h<b; i++){
        sum += h*f(a+i*h, p);
    }
    return sum;
}

double integrate_trapez_adap( double(*f)(double, void *), double a, double b, double acc, void *p){
    double sum_i, sum_ip1;  //zwei aufeinanderfolgende integralsummen
    double h;   //Schrittweite
    int N; //Anzahl der Schritte
    h = fabs(b-a)/4; //startschrittweite
    N = fabs(b-a)/h; //Startschrittzahl
    do{
        sum_i = sum_ip1;
        sum_ip1 = h/2*(f(a, p) + f(b, p));
        for(int i=1; i<N; i++){
            sum_ip1 += h*f(a+i*h, p);
        }

        h *= 0.5;
        N = fabs(b-a)/h;
        //printf("Der Momentane Wert des Integrals ist %15.6e \n", sum_ip1);
    }
    while(fabs(sum_i - sum_ip1) > acc);
    return sum_ip1;
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