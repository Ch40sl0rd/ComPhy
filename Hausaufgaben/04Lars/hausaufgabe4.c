//Lars Döpper
//make
// ./hausaufgabe4

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "arrayhelpers.h"
#include "filehelper.h"
#include "numerik_own.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double V_0; //parameter for unmodified potential
double v_0; //parameter for modified potential
double R;   //paramteter for unmodified potential
double mu; //parameter for exponential decay
/****************************************
 *  This function calculates the (m,n)-potential in the spacew room
 *  of the form:
 *  V(r) = V_0*((R/r)^m - m/n(R/r)^n)*n/(m-n)
 * 
 * \param r: value of r.
 * \param m: integer m
 * \param n: integer n
 * 
 * \return value of the (m,n) potential at r
 * **************************************/
double pot_space_norm(double r, int m, int n){
    return V_0*(n/(m-n))*( pow(R/r, m) - m/n*pow(R/r, n) );
}

/****************************************
 *  This function calculates the modified (m,n)-potential in the space room
 *  of the form:
 *  V(r) = v_0*((1/r)^m - m/n(1/r)^n)*n/(m-n)*exp(-mu*r)
 * 
 * \param r: value of r.
 * \param m: integer m
 * \param n: integer n
 * 
 * \return value of the (m,n) potential at r
 * **************************************/
double pot_space_modified(double r, int m, int n){
    if(mu < 0){
        printf("[pot_space_mod] wrong parameter for mu. Mu has be be greater than 0. Make Mu smaller than 0.\n");
        mu = -mu;
    }
    return v_0*n/(m-n)*(pow(1/r, m) - m/n*pow(1/r, n))*exp(-mu*r);
}

/****************************
 *  This function calculates the potiential in the
 *  impuls room for a potential of the form 
 *  V(r) = V_0*((R/r)^m - m/n(R/r)^n)*n/(m-n)
 *  in the space room by using numerical integration
 *  with the gaus-legendre-methode
 * 
 *  \param p: impuls p
 *  \param p_s: impuls p'
 *  \param m: integer m
 *  \param n: integer n
 * 
 *  \return value of the potential at p, p'
 * *************************/
double pot_impuls(double p, double p_s, int m, int n,
    double y_start, double y_end, int y_steps){

    double pot; //value of the potential;
    gsl_integration_glfixed_table *yw_table;
    size_t i;
    double *y_vals, *w_vals;
    double y, w, pot_space;
    y_vals = (double*)malloc(sizeof(double)*y_steps);
    w_vals = (double*)malloc(sizeof(double)*y_steps);

    //create table for gaus-legendre integration
    yw_table = gsl_integration_glfixed_table_alloc(y_steps);
    if(!yw_table){
        printf("[pot_impuls] Error while creating weights for gauss-legendre-intzegration.\n");
        exit(1);
    }
    //set points and weights for gauss-legendre integration
    for(i=0; i<y_steps; i++){
        gsl_integration_glfixed_point(y_start, y_end, i, &y_vals[i], &w_vals[i], yw_table);
    }
    
    pot = 0.0;
    for(i=0; i<y_steps;i++){
        y = y_vals[i];
        w = w_vals[i];
        pot_space = pot_space_modified(y, m, n);
        printf("y: %15.6e\t w:%15.6e\t pot: %15.6e\n", y, w, pot_space);
        pot += y*y*pot_space*sin(p*y)*sin(p_s*y)/(p*p_s*y*y);        
    }
    pot *= 2/M_PI;

    free(y_vals); free(w_vals);gsl_integration_glfixed_table_free(yw_table);
    return pot;
}

/********************************************************
 *  This function calculates the analytical value of the
 *  (2,1) potential in the impulse space for given p and p'
 * 
 *  \param p: value of p
 *  \param p_s: value of p'
 * 
 *  \return: analytical value of the (2,1) potential
 * *******************************************************/
double pot_impuls_analytical(double p, double p_s){
    return 0;
}

int main(int argc, char* argv[]){
    v_0 = 400.0;
    mu = 1;
    printf("Impuls bei p, p_s : %15.6e.\n", pot_impuls(10.0, 20.0, 2, 1, 0.0, 10.0, 100));
    return 0;
}