//Dirk Knott und Lars Döpper
//make
// ./hausaufgabe5

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "arrayhelpers.h"
#include "filehelper.h"

//Globale Variablen
int nt; //Anzahl der Zeitschritte
int nd; //Anzahl der Ortspunkte

int idx(int n, int j){
    if(j>=nd) j=nd; //Behandlung des rechten Randes bzgl außerhalb des Gebiets
    if(j<0) j=0; //Behandlung des rechten Rands in Orts-Richtung
    if(n>=nt || n<0){    //Behandlung der Zeit außerhalb des Gebietes, sollte nicht passieren
        n=nt; 
        printf("[idx] Das sollte nicht passieren.\n");
    }
    return nd*n + j;
}

void solve_kdv(double t_start, double t_end, double x_start, double x_end, double N, char* const filename){
    //Variablen deklarieren
    double d;   //Schrittweite in Zeit-Richtung
    double h;   //Schrittweite in Orts-Richtung
    double *u;  //Feld für die Werte der Funktion, ein 2d-array
    double *x;  //Feld für die x-Koordinaten
    double *t;  //Feld für die Zeit-Koordinaten
    double fak1, fak2;
    double *func_end;

    //Schrittweiten berechnen
    h = fabs((x_end-x_start)/(double)(nd-1));   
    d = fabs((t_end-t_start)/(double)(nt-1));

    //Konvergenzkriterium untersuchen
    if(d>h*h*h/2.6) printf("[solve_kdv] Falsche Schrittweite, Konvergenz nicht garantiert.\n");

    //Speicher für u allokieren
    u = (double*)malloc(sizeof(double)*nt*nd);
    x = (double*)malloc(sizeof(double)*nd);
    t = (double*)malloc(sizeof(double)*nt);
    func_end = (double*)malloc(sizeof(double)*nd);

    //Startwerte initialisieren:
    x[0] = x_start;
    u[idx(0,0)] = 0.0;
    x[nd-1] = x_end;
    u[idx(0, nd-1)] = 0.0;
    for(int i=1; i<nd-1; i++){
        x[i] = x_start+i*h;
        
        u[idx(0,i)]= (double)(-N*(N+1))/pow(cosh(x[i]), 2);
    }
    t[0] = t_start;

    //ersten Zeitschritt mit anderen Differenzenquotienten durchführen.
    t[1] = t_start + d;
    for(int j=0; j<nd; j++){
        fak1 = 6.0*u[idx(0, j)]*(u[idx(0,j+1)]-u[idx(0,j-1)])/(2.0*h);
        fak2 = -(u[idx(0,j+2)] - 2.0*u[idx(0,j+2)] + 2.0*u[idx(0,j-1)] - u[idx(0,j-2)])/(2.0*h*h*h);
        u[idx(1,j)] = d*(fak1+fak2) + u[idx(0,j)];
    }
    

    for(int n=1; n<nt-1; n++){
        t[n] = t_start + n*d;
        for(int j=0; j<nd; j++){
            //Brechnung der neuen Punktes u_n^j
            //Setze den Rand immer auf 0
            if(j==0 || j==nd-1){
                u[idx(n+1,j)] = 0.0;
            }
            else{
                fak1 = ((u[idx(n,j+1)] + u[idx(n,j)] + u[idx(n,j-1)])/3.0) * ((u[idx(n, j+1)] - u[idx(n,j-1)])/(2.0*h));
                fak2 = (u[idx(n,j+2)] - 2.0*u[idx(n,j+1)] + 2.0*u[idx(n,j-1)] - u[idx(n,j-2)])/(2.0*h*h*h);
                u[idx(n+1, j)] = 2.0*d*(6.0*fak1 - fak2) + u[idx(n-1,j)];
                
            }
        }
    }

    for(int i=0; i<nd; i++){
        func_end[i] = u[idx(nt-1, i)];
    }

    print_data_table_double(filename, x, func_end, nd);
    //print_data2file(filename2, x, nd);
    //print_data2file(filename3, t, nt);
    free(u); free(x); free(t);free(func_end);
}

int main(){
    double x_start, x_end, t_start, t_end;
    nt = 25000;
    nd = 1201;
    x_start = -30.0;
    x_end = 30.0;
    t_start = 0.0; 
    t_end = 1.0;
    solve_kdv(t_start, t_end, x_start, x_end, 2, "funktion2.txt");
    solve_kdv(t_start, t_end, x_start, x_end, 1, "funktion1.txt");
    nd=1001;
    nt = 25000;
    x_start = -5.0;
    x_end = 45.0;
    solve_kdv(t_start, t_end, x_start, x_end, 3, "funktion3.txt");
}