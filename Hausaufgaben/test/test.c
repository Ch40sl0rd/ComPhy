#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "module1.h"
#include "arrayhelpers.h"
#include "numerik_own.h"

double function(double x, void* p){
    return sin(x);
}

int main(int argc, char** argv){
    double a,b;
    int n;
    if(argc < 4){
        printf("Verwende standart werte. \n");
        a = 0.;
        b = 5.;
        n = 9;
    }
    else{
        a = atof(argv[1]);
        b = atof(argv[2]);
        n = atoi(argv[3]);
    }
    printf("Das ist ein test.\n");
    //print1();
    //print2();
    
    double* data = create_array_double(a, b, n);
    double* y_data = create_array_function_double(data, function, n, NULL);
    array_ausgabe_double(y_data, n, 1);
    double *deriv = derivate_sym_one_array(y_data, a, b, n);
    array_ausgabe_double(deriv, n, 1);
    double* sec_derv = derivate_sym_two_array(y_data, a, b, n);
    array_ausgabe_double(sec_derv, n, 1);
    
}
