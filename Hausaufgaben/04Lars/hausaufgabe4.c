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

int dim; //dimension of the matrix, nummer der stützstellen für p
int maxiter=6; // maximum number of iteration steps

int m;  //parameter for potential
int n;  //parameter for potential

double y_start; //start value for integration of potential
double y_end;   //end value for integration of potential
int y_steps;    //number of values for y for integration

double p_start;
double p_end;


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
double pot_space_norm(double r){
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
double pot_space_modified(double r){
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
 *  v(y) = v_0*((1/y)^m - m/n(1/y)^n)*n/(m-n)
 *  in the space room by using numerical integration
 *  with the gaus-legendre-methode
 * 
 *  \param p: impuls p
 *  \param p_s: impuls p'
 *  
 *  \return value of the potential at p, p'
 * *************************/
double pot_impuls(double p, double p_s){
    
    if(p==0.0){
        p=1e-5;
    }
    if(p_s==0.0){
        p_s=1e-5;
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
        pot_space = pot_space_modified(y);
        //printf("y: %15.6e\t w:%15.6e\t pot: %15.6e\n", y, w, pot_space);
        pot += w*y*y*pot_space*sin(p*y)*sin(p_s*y)/(p*p_s*y*y);        
    }
    pot *= 2/M_PI;

    free(y_vals); free(w_vals);gsl_integration_glfixed_table_free(yw_table);
    return pot;
}

/**
 * Diese Methode erstellt eine Liste der Potentiale im Ortsraum mit
 * p_start <= p <= p_max mit dim-vielen Schritten für ein festes p'.
 */
void pot_impuls_list(double p_s){
    double pot; //value of the potential;
    double p;
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
    
   
    for(int j=0; j<dim; j++){
        pot = 0.0;
        p = p_array[j];
        for(i=0; i<y_steps;i++){
            y = y_vals[i];
            w = w_vals[i];
            pot_space = pot_space_modified(y);
            //printf("y: %15.6e\t w:%15.6e\t pot: %15.6e\n", y, w, pot_space);
            pot += w*y*y*pot_space*sin(p*y)*sin(p_s*y)/(p*p_s*y*y);        
        }
        pot *= 2/M_PI;
        V_array[j] = pot;
    }
    free(y_vals); free(w_vals); gsl_integration_glfixed_table_free(yw_table);
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
    if(p==0.0) p=1e-5;
    if(p_s==0.0) p_s=1e-5;
    double pre_fak, fak1, fak2, fak3;
    pre_fak = v_0/(2.0*M_PI*p*p_s);
    fak1 = 2.0*(p_s-p)*atan((p-p_s)/mu);
    fak2 = 2.0*(p_s+p)*atan((p+p_s)/mu);
    fak3 = (2.0+mu)*(log(1.0+(p-p_s)*(p-p_s)/(mu*mu)) - log(1.0+ (p_s+p)*(p_s+p)/(mu*mu)));
    return pre_fak*(fak1+fak2+fak3);
}

void prep_trapez(double *p, double *w){
    double h = (p_end-p_start)/(dim-1);
    for(int i=0; i<dim; i++){
        p[i] = p_start +i*h;
        w[i] = h;
    }
    w[0] = h/2.0;
    w[dim-1] = h/2.0;
}

void prep_mat(){
    //allokiere speicher für die stützstellen von p, w und V
    p_array = (double*)malloc(sizeof(double)*dim); //array für impulspunkte
    w_array = (double*)malloc(sizeof(double)*dim); //array für gewichte der trapez-integration
    V_array = (double*)malloc(sizeof(double)*dim*dim);  // matrix für werte des Potentials

    prep_trapez(p_array, w_array);

    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            V_array[i + j*dim] = pot_impuls(p_array[i], p_array[j])*p_array[j]*p_array[j]*w_array[j];
        }
    }
}

void prep_mat_neu(){
    double pi, pj;
    double pot;
    //allokiere speicher für die stützstellen von p, w und V
    p_array = (double*)malloc(sizeof(double)*dim); //array für impulspunkte
    w_array = (double*)malloc(sizeof(double)*dim); //array für gewichte der trapez-integration
    V_array = (double*)malloc(sizeof(double)*dim*dim);  // matrix für werte des Potentials

    prep_trapez(p_array, w_array);

    //Berechne das Potential für alle P aus der Liste p_array
    //Mit gaus-legendre-integration. Dafür benötigen wir die Hilfsfelder
    //ingesamt nur ein mal, da sich an den Grenzen für y nichts ändert.
    //Potential für pi, pj;
    
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
    
    
    //iteriere über alle Matrix elemente:
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            pot = 0.0;
            pi = p_array[i];
            pj = p_array[j];
            for(int k=0; k<y_steps;k++){
                y = y_vals[k];
                w = w_vals[k];
                pot_space = pot_space_modified(y);
                //printf("y: %15.6e\t w:%15.6e\t pot: %15.6e\n", y, w, pot_space);
                pot += w*y*y*pot_space*sin(pi*y)*sin(pj*y)/(pi*pj*y*y);        
            }
            pot *= 2/M_PI; //fertiges potential im Impulsraum für pi, pj
            V_array[i +j*dim] = pot*pj*pj*w_array[j]; //wert des Potentials mit Vorfaktoren für Arnoldi-Verfahren.
        }
    }
    free(y_vals); free(w_vals);gsl_integration_glfixed_table_free(yw_table);
}

void apply_mat(double *vec_in, double *vec_out, double E){
    double E_fakt; //Vorfaktor für Matrixelemente
    for(int i=0; i<dim; i++){
        vec_out[i] = 0.0;
        E_fakt = 1.0/(E-0.5*p_array[i]*p_array[i]);
        for(int j=0; j<dim; j++){
            vec_out[i] += E_fakt*V_array[i*dim +j]*vec_in[j]; //Die restliche Faktoren liegen schon im Pot-Array
        }
    }
}

double eigenval(double E){
    double maxlambda, maxlambda_last; //Variablen für die letzen beiden Eigenwerte
    //wird als Konvergenzkriterium verwendet
    double  *a_mat; //Matrix A
    double *v_vecs;  //Feld  der Basisvektoren v_i
    double *w_vec;  //Feld für Hilfsvektoren w_{j+1}
    double *w_tilde_vec; //Feld  für normierte Hilfsvektoren w_tilde_{j+1}
    double *v_help;  //Hilfsvektoren für Zwischenergebnis
    double *VL, *VR;   //Speicher für lapack-bib
    double *c;  //Feld für die Eigenvektoren
    double *WR, *WI;    //Feld für Real- und Imanginärteil der Eigenwerte
    double *work;   //Hilsfeld  für  lapack-bib
    double norm, dot_prod_val;  //Hilfswerte für Norm und Skalarprodukt
    char jobvl='N',jobvr='V'; 

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

    //initianlisiere matrix a_mat als 0
    for(int i=0; i<dim_red*dim_red; i++){
        a_mat[i] = 0.0;
    }

    //setze den startvektor als konstant 1 fest.
    for(int i=0; i<dim; i++){
        v_vecs[i + 0*dim_red] = 1.0;
    }
    maxlambda = 0.0;
    norm = sqrt(dot_prod(&v_vecs[0*dim], &v_vecs[0*dim]));
    //normalisiere den neuen startvektor
    multiply_vec(1/norm, &v_vecs[0*dim], &v_vecs[0*dim]);

    for(int k =0; k<maxiter; k++){
        apply_mat(&v_vecs[k*dim], w_vec, E);
        //Gram-Schmidt-Orthogonalisierung######################################
        //übertrage w_vec in w_tilde_vec
        multiply_vec(1.0, w_vec, w_tilde_vec);
        for(int j=0; j<=k; j++){
            dot_prod_val = dot_prod(&v_vecs[j*dim], w_vec);
            a_mat[j+dim_red*k] = dot_prod_val;
            multiply_vec(-dot_prod_val, &v_vecs[j*dim], v_help );
            addvec(w_tilde_vec, v_help, w_tilde_vec);
        }
        a_mat[k+1+dim_red*k] = sqrt(dot_prod(w_tilde_vec, w_tilde_vec));
        norm = 1.0/a_mat[k+1+dim_red*k];
        multiply_vec(norm, w_tilde_vec, &v_vecs[(k+1)*dim]);
    }

    //Wir haben jetzt die Matrix A bestimmt und diagonalisieren diese jetzt.
    dim--;
    dgeev_(&jobvl, &jobvr, &dim_red, a_mat, &lda, WR, WI, VL, &ldvl, VR, &ldvr, work, &lwork, &info);

    maxlambda = 0.0;
    for(int i=0; i<dim_red; i++){
        if(fabs(WI[i]<1e-6)){
            maxlambda = fmax(maxlambda, WR[i]);
        }
    }
    if(info!=0){
        printf("Errorcode: %d\n", info);
        abort();
    }

    free(a_mat); free(v_vecs); free(w_vec); free(w_tilde_vec); free(v_help); free(VL); free(VR);
    free(c); free(WR); free(WI); free(work);
    return maxlambda;
}

double function_search(double E){
    return eigenval(E) - 1.0;
}

int main(int argc, char* argv[]){
    //###################### Aufgabe 9.3 ######################################
    double *V_error, *V_analyt, *V_error_max, *n_y_array;
    p_start = 0.0; p_end = 200.0;
    dim = 400;
    y_start = 0.0; y_end = 10.0; y_steps=1000;
    v_0 = 400.0;
    mu = 1.0;
    R = 1.0;
    double p_tilde = 1.0;
    m = 2; n=1;
    p_array = create_array_double(p_start, p_end, dim);
    V_analyt = (double*)malloc(sizeof(double)*dim);
    V_error = (double*)malloc(sizeof(double)*dim);
    V_array = (double*)malloc(sizeof(double)*dim);
    n_y_array = (double*)malloc(sizeof(double)*4);
    V_error_max = (double*)malloc(sizeof(double)*4);

    for(int i=2; i<6; i++){
        y_steps = pow(10.0, i);
        n_y_array[i-2] = y_steps;
        pot_impuls_list(p_tilde);
        double max_error = 0.0;
        for(int j=0; j<dim; j++){
            V_analyt[j] = pot_impuls_analytical(p_array[j], p_tilde);
            V_error[j] = fabs((V_array[j]-V_analyt[j])/V_analyt[j]);
            max_error = fmax(max_error, V_error[j]);
        }
        printf("Maximaler Fehler für ny=%d ist %15.6e\n", y_steps, max_error);
        V_error_max[i-2] = max_error;
    }
    print_data_table_double(NULL, n_y_array, V_error_max, 4);
    free(p_array); free(V_array); free(V_analyt); free(V_error); free(n_y_array); free(V_error_max);

    //###################### Aufgabe 9.6 ########################################
    //Wir ändern ein paar parameter für das Lennard-Jones Potential.
    m = 12; n=6;
    p_end = 400.0; dim = 2000;
    y_start = 0.4; y_end=10.0; y_steps = 20000;
    p_array = create_array_double(0.0, p_end, dim);
    V_array = (double*)malloc(sizeof(double)*dim);
    pot_impuls_list(p_tilde);
    print_data_table_double("Lennard_jones_pot.txt", p_array, V_array, dim);
    free(p_array); free(V_array);

    
    return 0;
}