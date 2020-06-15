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
#define FINT int

double V_0; //parameter for unmodified potential
double v_0; //parameter for modified potential
double R;   //paramteter for unmodified potential
double mu; //parameter for exponential decay

double *p_array;
double *w_array;
double *V_array;

int dim; //dimension of the matrix
int maxiter=5; // maximum number of iteration steps

/* dgeev_ als externe Funktion deklarieren 
   Die hier gegebenen Parameter stimmen mit der Definition
   in LAPACK Bibliothek ueberein. 
   Uebergabe der Adresse eines Objekts auch fuer nicht-Felder 
    (in FORTRAN ist das Standard)       
   FORTRAN Name DGEEV    ->    dgeev_                      */

extern void dgeev_(char *jobvl,char *jobvr,FINT *n,double *a,FINT *lda,
                   double *wr,double *wi,double *vl,FINT *ldvl, 
                   double *vr,FINT *ldvr, 
                   double *work,FINT *lwork,FINT *info);

/* Routine fuer Nullstellensuche mit dem Sekantenverfahren  (siehe Beispiel 3.3-secant.c) */ 
double secant(double x1, double x2, double (*func)(double), int *schritt)
/* x1,x2     Startwerte 
   func      ist die "Referenz" auf eine Funktion mit einem double Parameter, 
             die double zurueckgibt (Referenz = Adresse der Funktion)  
   schritt   ist auch Referenz: Veraenderungen an der Variable wirken sich auf 
                                das aufrufende Programm aus !!! */
{
  const double tol=1e-12; /* geforderte Genauigkeit als Konstante */
  double xn;              /* neuer Schaetzwert */ 

  *schritt=0;  /* noch kein Schritt */

  do
    {         /* naechster Schaetzwert x1,x2 -> xn*/
      xn=x2-func(x2)*(x2-x1)/(func(x2)-func(x1));
      x1=x2;     /* bereite den naechsten Schritt vor:  x2 -> x1 */  
      x2=xn;     /*                                     xn -> x2 */

      (*schritt)++;      /* Schritte=Schritte+1 */
    }
  while(fabs(x2-x1)>tol);   /* solange Genauigkeitsziel nicht erreicht */

  return xn;   /* Gebe Nullstelle zurueck */
}

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

void  addvec(double *vin1, double *vin2, double *vout){
    for(int i=0; i<dim; i++){
        vout[i] = vin1[i] + vin2[i];
    }
}

void multiply_vec(double a, double *vin, double *vout){
    for(int i=0; i<dim; i++){
        vout[i] = a*vin[i];
    }
}

double dot_prod(double *vin1, double *vin2){
    double res =0.0;
    for(int i=0; i<dim; i++){
        res += vin1[i]*vin2[i];
    }
    return res;
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
    
    if(p==0.0 || p_s==0.0){
        return 1e-5;
    }

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

double eigenval(double E){
    double maxlambda;
    double  *a_mat; //Matrix A
    double *v_vecs;  //Feld  der Basisvektoren v_i
    double *w_vec;  //Feld für Hilfsvektoren w_{j+1}
    double *w_tilde_vec; //Feld  für normierte Hilfsvektoren w_tilde_{j+1}
    double *v_help;  //Hilfsvektoren für Zwischenergebnis
    double *VL, *VR;   //Speicher für lapack-bib
    double *c;  //Feld für die Eigenvektoren
    double *WR, *WI;    //Feld für Real- und Imanginärteil der Eigenwerte
    double *work;   //Hilsfeld  für  lapack-bib
    double norm, dot_prod;  //Hilfswerte für Norm und Skalarprodukt

    //Hilfsvariablen  für die Verwendung von Lapack mit Fortran
    FINT info;
    FINT lda = maxiter+1, ldvl = maxiter+1, ldvr=maxiter+1;
    FINT lwork = 4*maxiter;
    FINT dim_red = maxiter+1;

    //Speicherallokierung für die verschiedenen Felder mit reduzierter Dimension
    a_mat = (double*)malloc(sizeof(double)*dim_red*dim_red);
    v_vecs = (double*)malloc(sizeof(double)*dim_red*dim);
    w_vec = (double*)malloc(sizeof(double)*dim);
    w_tilde_vec = (double*)malloc(sizeof(double)*dim);
    v_help = (double*)malloc(sizeof(double)*dim);
    WR = (double*)malloc(sizeof(double)*dim_red);
    WI = (double*)malloc(sizeof(double)*dim_red);
    work = (double*)malloc(sizeof(double)*lwork);
    c = (double*)malloc(sizeof(double)*dim_red);
    VR = (double*)malloc(sizeof(double)*dim_red*dim_red);
    VL = (double*)malloc(sizeof(double)*dim_red*dim_red);

}

double function_search(double E){
    return eigenval(E) - 1.0;
}

int main(int argc, char* argv[]){
    v_0 = 400.0;
    mu = 1;
    printf("Impuls bei p, p_s : %15.6e.\n", pot_impuls(10.0, 20.0, 2, 1, 0.0, 10.0, 100));
    return 0;
}