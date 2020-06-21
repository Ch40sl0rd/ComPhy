//Lars Döpper
//make
// ./hausaufgabe4 [variation der Parameter, 1 für ja, 0 für nein]
//Es empfiehlt sich, die Variation der Parameter für den ersten Durchlauf auszukommentieren,
//da das Programm auf meinem Rechner gut eineinhalb Stunden läuft.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "arrayhelpers.h"
#include "filehelper.h"
#include<ctype.h>
#include<sys/time.h>

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
    //Wir prüfen einen divide by zero error
    if(r==0.0) r=pow(10.0, -(7-(double)m/2.0));
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
    if(p_s==0.0) p_s=1e-5;
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
        if(p==0.0) p=1e-5; //check for divide by zero
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

void prep_mat_anal(){
    double pi, pj;
    //allokiere speicher für die stützstellen von p, w und V
    p_array = (double*)malloc(sizeof(double)*dim); //array für impulspunkte
    w_array = (double*)malloc(sizeof(double)*dim); //array für gewichte der trapez-integration
    V_array = (double*)malloc(sizeof(double)*dim*dim);  // matrix für werte des Potentials

    prep_trapez(p_array, w_array);

    for(int i=0; i<dim; i++){
        pi = p_array[i];
        if(pi==0.0) pi=1e-5;
        for(int j=0; j<dim; j++){
            pj = p_array[j];
            if(pj==0.0) pj=1e-5;
            V_array[i + j*dim] = pot_impuls_analytical(pi, pj)*pj*pj*w_array[j];
        }
    }
}

/**
 *  Diese Methode bereit die Potentialmatrix für ein m,n Potential vor
 *  Dabei wird der speicher für die listen p, w und V neu allokiert.
 * 
 */
void prep_mat_neu(){
    double pi, pj;
    double pot;
    //allokiere speicher für die stützstellen von p, w und V
    p_array = (double*)malloc(sizeof(double)*dim); //array für impulspunkte
    w_array = (double*)malloc(sizeof(double)*dim); //array für gewichte der trapez-integration
    V_array = (double*)malloc(sizeof(double)*dim*dim);  // matrix für werte des Potentials

    prep_trapez(p_array, w_array);  //vorbereitung der Gewichte

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
                if(pi == 0.0) pi = 1e-5;
                if(pj == 0.0) pj = 1e-5;
                y = y_vals[k];
                w = w_vals[k];
                pot_space = pot_space_modified(y);
                //printf("y: %15.6e\t w:%15.6e\t pot: %15.6e\n", y, w, pot_space);
                pot += w*pot_space*sin(pi*y)*sin(pj*y)/(pi*pj);        
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
        if(E_fakt==0.0){
            E_fakt = 1e-6;
        }
        for(int j=0; j<dim; j++){
            vec_out[i] += E_fakt*V_array[i*dim +j]*vec_in[j]; //Die restliche Faktoren liegen schon im Pot-Array
        }
    }
}

double eigenval(double E){
    double maxlambda; //Variablen für die letzen beiden Eigenwerte
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
    dim_red--;
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

    double timeused;                /* Variabeln zur Zeitmessung */
    struct timeval tvstart,tvend;
    double *free_param, *rel_error;
    int steps_param;
    int verbose;
    if(argc>=2){
        verbose = atoi(argv[1]);
    }
    else verbose =0;
    if(verbose){
        printf("[main] Es wird eine Variation der Parameter durchgeführt, dies kann sehr lange dauern.\n");
    }
    else{
        printf("[main] Es wird keine Variation der Parameter durchgeführt.\n");
    }
    //printf("verbose: %d\n", verbose);
    gettimeofday(&tvstart,NULL);
    //###################### Aufgabe 9.3 ######################################
    printf("[main] Wir berechnen jetzt numerisch das (2,1) Potential im Impulsraum\nund vergleichen dieses mit der analytischen Lösung.\n");
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

    for(int i=2; i<5; i++){ //hier noch wieder 6 als Obergrenze einfügen.
        y_steps = pow(10.0, i);
        n_y_array[i-2] = y_steps;
        pot_impuls_list(p_tilde);
        double max_error = 0.0;
        for(int j=0; j<dim; j++){
            V_analyt[j] = pot_impuls_analytical(p_array[j], p_tilde);
            V_error[j] = fabs((V_array[j]-V_analyt[j])/V_analyt[j]);
            max_error = fmax(max_error, V_error[j]);
        }
        //printf("[main] Maximaler Fehler für ny=%d ist %15.6e\n", y_steps, max_error);
        V_error_max[i-2] = max_error;
    }
    print_data_table_double("9_3_variation_stützstellen.txt", n_y_array, V_error_max, 4);
    free(p_array); free(V_array); free(V_analyt); free(V_error); free(n_y_array); free(V_error_max);

    //###################### Aufgabe 9.4 ########################################
    //Zunächst ein paar Konstanten
    printf("[main] Aufgabe 9.4 für die Energieeigenwerte.\n");
    y_steps = 1000;
    maxiter = 50;
    prep_mat_neu();
    //printf("[main] Matrix ist fertig bestimmt.\n");
    int steps_num = 11;
    double energy_num, energy_anal;
    energy_num = secant(-100.0, -10.0, function_search, &steps_num);
    printf("Der Energieeigenwert numerisch lautet: %15.6e nach %d Schritten\n", energy_num, steps_num);
    free(p_array); free(V_array); free(w_array);
    prep_mat_anal();
    //printf("Die analytische Matrix ist fertig bestimmt.\n");
    energy_anal = secant(-100.0, -10.0, function_search, &steps_num);
    printf("Der Energieeigenwert analytisch lautet: %15.6e nach %d Schritten\n", energy_anal, steps_num);
    free(p_array); free(w_array); free(V_array);
    
    if(verbose){
    //Wir variieren jetzt einige Paramter und betrachten die Veränderung der Eigenenergie.
    
    steps_param=7;
    //zunächst varrieren wir die Anzahl der y-Stützstellen
    free_param = (double*)malloc(sizeof(double)*steps_param);
    rel_error = (double*)malloc(sizeof(double)*steps_param);
    printf("[main] Variation der y-Stützstellen\n");
    y_steps = 100;
    for(int i=0; i<steps_param; i++){
        y_steps *= 2;
        free_param[i] = y_steps;
        prep_mat_neu();
        energy_num = secant(-100.0, -10.0, function_search, &steps_num);
        //printf("Der Energieeigenwert numerisch lautet: %15.6e nach %d Schritten\n", energy_num, steps_num);
        free(p_array); free(V_array); free(w_array);
        prep_mat_anal();
        energy_anal = secant(-100.0, -10.0, function_search, &steps_num);
        //printf("Der Energieeigenwert analytisch lautet: %15.6e nach %d Schritten\n", energy_anal, steps_num);
        rel_error[i] = fabs((energy_anal-energy_num)/energy_anal);
        free(p_array); free(w_array); free(V_array);
    }
    print_data_table_double("variation_ny.txt", free_param, rel_error, steps_param);
    free(rel_error); free(free_param);

    //Jetzt varrieren wir den Maximimalen Impuls bei ny = 1000
    y_steps = 1000;
    steps_param = 10;
    free_param = (double*)malloc(sizeof(double)*steps_param);
    rel_error = (double*)malloc(sizeof(double)*steps_param);
    printf("[main] Variation des maximalen Impulswertes.\n");
    for(int i=0; i<steps_param; i++){
        p_end = 200.0 +i*20.0;
        free_param[i] = p_end;
        prep_mat_neu();
        energy_num = secant(-100.0, -10.0, function_search, &steps_num);
        //printf("Der Energieeigenwert numerisch lautet: %15.6e nach %d Schritten\n", energy_num, steps_num);
        free(p_array); free(V_array); free(w_array);
        prep_mat_anal();
        energy_anal = secant(-100.0, -10.0, function_search, &steps_num);
        //printf("Der Energieeigenwert analytisch lautet: %15.6e nach %d Schritten\n", energy_anal, steps_num);
        rel_error[i] = fabs((energy_anal-energy_num)/energy_anal);
        free(p_array); free(w_array); free(V_array);
    }
    print_data_table_double("variation_pmax.txt", free_param, rel_error, steps_param);
    free(rel_error); free(free_param);

    //Als letztes varriieren wir die Anzahl der Impuls-Stützstellen np.
    p_end = 200.0;
    steps_param = 10;
    free_param = (double*)malloc(sizeof(double)*steps_param);
    rel_error = (double*)malloc(sizeof(double)*steps_param);
    printf("[main] Variation der p-Stützstellen\n");
    for(int i=0; i<steps_param; i++){
        dim = 400 +i*20;
        free_param[i] = dim;
        prep_mat_neu();
        energy_num = secant(-100.0, -10.0, function_search, &steps_num);
        //printf("Der Energieeigenwert numerisch lautet: %15.6e nach %d Schritten\n", energy_num, steps_num);
        free(p_array); free(V_array); free(w_array);
        prep_mat_anal();
        energy_anal = secant(-100.0, -10.0, function_search, &steps_num);
        //printf("Der Energieeigenwert analytisch lautet: %15.6e nach %d Schritten\n", energy_anal, steps_num);
        rel_error[i] = fabs((energy_anal-energy_num)/energy_anal);
        free(p_array); free(w_array); free(V_array);
    }
    print_data_table_double("variation_np.txt", free_param, rel_error, steps_param);
    free(rel_error); free(free_param);
    }

    

    


    //###################### Aufgabe 9.6 ########################################
    //Wir ändern ein paar parameter für das Lennard-Jones Potential.
    printf("[main] Berechnung des Lennard-Jones-Potentials im Impulsraum.\n");
    m = 12; n=6; mu=0.0;
    p_end = 400.0; dim = 2000;
    y_start = 0.4; y_end=10.0; y_steps = 20000;
    p_array = create_array_double(0.0, p_end, dim);
    V_array = (double*)malloc(sizeof(double)*dim);
    pot_impuls_list(p_tilde);
    print_data_table_double("Lennard_jones_pot.txt", p_array, V_array, dim);
    free(p_array); free(V_array);

    //##################### Aufgabe 9.7 #########################################
    //Berechnung des Energy-Eigenwertes für Lennard-Jones Potential.
    printf("[main] Berechnung des Energieeigenwerts für das Lennard-Jones-Potential.\n");
    v_0 = 400.0; R=1.0; p_start = 0.0;
    y_steps = 5000; dim = 500; p_end=400.0;
    y_start = 0.4; y_end = 10.0;
    m = 12; n=6;
    prep_mat_neu();
    //printf("[main] Die Matrix ist fertig bestimmt.\n");
    energy_num = secant(-100.0, -50.0, function_search, &steps_num);
    printf("Der Energieigenwert lautet %15.6e nach %d Schritten.\n", energy_num, steps_num);
    free(p_array); free(V_array); free(w_array);
    
    if(verbose){
    //Wir variieren jetzt einige parameter und beobachten den Eigenwert.
    steps_param = 10;
    //Zunächst den minimalwert für y
    printf("[main] Variation des Minimalwerts y_min.\n");
    y_steps = 5000; dim = 500; p_end = 400.0;
    free_param = (double*)malloc(sizeof(double)*steps_param);
    rel_error = (double*)malloc(sizeof(double)*steps_param);
    for(int i=0; i<steps_param; i++){
        y_start = 0.4 +i*0.05;
        free_param[i] = y_start;
        prep_mat_neu();
        energy_num = secant(-100.0, -50.0, function_search, &steps_num);
        //printf("Der Energieigenwert lautet %15.6e nach %d Schritten.\n", energy_num, steps_num);
        rel_error[i] = energy_num;
        free(p_array); free(V_array); free(w_array);
    }
    print_data_table_double("lj_variation_ymin.txt", free_param, rel_error, steps_param);
    free(rel_error); free(free_param);

    //Variation des maximalen Impulses.
    printf("[main] Variation des Maximalen Impulses p_max.\n");
    y_start =0.4;
    y_steps = 5000; dim = 500; p_end = 400.0;
    free_param = (double*)malloc(sizeof(double)*steps_param);
    rel_error = (double*)malloc(sizeof(double)*steps_param);
    for(int i=0; i<steps_param; i++){
        p_end = 400 +i*60;
        free_param[i] = p_end;
        prep_mat_neu();
        energy_num = secant(-100.0, -50.0, function_search, &steps_num);
        //printf("Der Energieigenwert lautet %15.6e nach %d Schritten.\n", energy_num, steps_num);
        rel_error[i] = energy_num;
        free(p_array); free(V_array); free(w_array);
    }
    print_data_table_double("lj_variation_pmax.txt", free_param, rel_error, steps_param);
    free(rel_error); free(free_param);

    //Variation von Anzahl Impulsstellen
    printf("[main] Variation der Stützstellen np.\n");
    y_steps = 5000; dim = 200; p_end = 400.0; y_start = 0.4;
    free_param = (double*)malloc(sizeof(double)*steps_param);
    rel_error = (double*)malloc(sizeof(double)*steps_param);
    for(int i=0; i<steps_param; i++){
        dim = 200 +i*100;
        free_param[i] = dim;
        prep_mat_neu();
        energy_num = secant(-100.0, -50.0, function_search, &steps_num);
        //printf("Der Energieigenwert lautet %15.6e nach %d Schritten.\n", energy_num, steps_num);
        rel_error[i] = energy_num;
        free(p_array); free(V_array); free(w_array);
    }
    print_data_table_double("lj_variation_np.txt", free_param, rel_error, steps_param);
    free(rel_error); free(free_param);

    //Variation von Anzahl Integrationsstellen
    printf("[main] Variation der Stützstellen der Integration ny.\n");
    y_steps = 1000; dim = 200; p_end = 400.0; y_start = 0.4;
    free_param = (double*)malloc(sizeof(double)*steps_param);
    rel_error = (double*)malloc(sizeof(double)*steps_param);
    for(int i=0; i<steps_param; i++){
        y_steps = 1000 + i*1000;
        free_param[i] = y_steps;
        prep_mat_neu();
        energy_num = secant(-100.0, -50.0, function_search, &steps_num);
        //printf("Der Energieigenwert lautet %15.6e nach %d Schritten.\n", energy_num, steps_num);
        rel_error[i] = energy_num;
        free(p_array); free(V_array); free(w_array);
    }
    print_data_table_double("lj_variation_ny.txt", free_param, rel_error, steps_param);
    free(rel_error); free(free_param);
    }

    /* Wie lange hat das in sec gedauert ? */
    gettimeofday(&tvend,NULL);
    timeused=tvend.tv_sec-tvstart.tv_sec;
    timeused=timeused+(tvend.tv_usec-tvstart.tv_usec)*1e-6;  /* Zeitdifferenz  in sec  */

    printf("[main] Das Programm hat %14.6le Sekunden benötigt.\n", timeused);


    
    return 0;
}