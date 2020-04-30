#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_integration.h>
#include <math.h>

void gaus_legendre(double a, double b, double *x_vals, double *w_vals, size_t n){
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

double gl_integrate(double a, double b, double *x, double *w, size_t n, double(*f)(double, void *), void *p){
    double sum;
    gaus_legendre(a,b, x,w, n);
    sum = 0.;
    for(int i=0; i<n; i++){
        sum += f(x[i], p)*w[i];
    }
    return sum;
}

double f(double x, void *p){
    return cos(log(x)*x)/x;
}

int main(){
    double a,b;
    size_t n;
    double exact, diff, sum;
    double *p;
    a = 1.;
    b = 500000;
    n = 500;

    p = (double*)malloc(sizeof(double)*2);
    p[0] = 2.;
    p[1] = 1.;

    double *x_vals = (double*)malloc(sizeof(double)*n);
    double *w_vals = (double*)malloc(sizeof(double)*n);

    
    sum = gl_integrate(a, b, x_vals, w_vals, n, &f, (void*)p);
    printf("Das Ergebnis lautet %15.6e\n", sum);

    free(x_vals);
    free(w_vals);
}