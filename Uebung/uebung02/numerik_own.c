#include <stdio.h>
double derivate_sym_one(double x, double h, double(*func)(double, void*), void* p){
    return (func(x+h, p)-func(x-h, p))/(2*h);
}

double derivate_sym_two(double x, double h, double(*func)(double, void*), void* p){
    return (func(x+h, p) + func(x-h, p) -2*func(x, p))/(h*h);
}
double derivate_sym_three(double x, double h, double(*func)(double, void*), void*p){
    return ( (func(x+2*h, p) - func(x-2*h, p)) -2*(func(x+h, p) - func(x-h, p)) )/(2*h*h*h);
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