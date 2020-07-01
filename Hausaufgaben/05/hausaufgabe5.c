//Dirk Knott und Lars Döpper
//make
// ./hausaufgabe5

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "arrayhelpers.h"
#include "filehelper.h"

int idx(int n, long nt, int j, long nd){
    if(j>=nd) j=nd; //Behandlung des rechten Randes bzgl außerhalb des Gebiets
    if(j<0) j=0; //Behandlung des rechten Rands in Orts-Richtung
    if(n>=nt || n<0){    //Behandlung der Zeit außerhalb des Gebietes, sollte nicht passieren
        n=nt; 
        printf("[idx] Das sollte nicht passieren.\n");
    }
    return nd*n + j;
}

void solve_kdv(double t_start, double t_end, double d, double x_start, double x_end, double h, double N, char* const filename1, char* const filename2){
    //Variablen deklarieren
    long nt;
    long nd;
    double *u;  //Feld für die Werte der Funktion, ein 2d-array
    double *x;  //Feld für die x-Koordinaten
    double *t;  //Feld für die Zeit-Koordinaten
    double fak1, fak2;
    double *func_end;

    nt = (long)(fabs(t_end-t_start)/d +1);
    nd = (long)(fabs(x_end-x_start)/h +1);
    
    //printf("[solve_kdv] Parameter d=%lf und nd=%lf\n", d, h);

    //Konvergenzkriterium untersuchen
    if(d>h*h*h/2.6) printf("[solve_kdv] Falsche Schrittweite, Konvergenz nicht garantiert.\n");

    //Speicher für u allokieren
    u = (double*)malloc(sizeof(double)*nt*nd);
    x = (double*)malloc(sizeof(double)*nd);
    t = (double*)malloc(sizeof(double)*nt);
    func_end = (double*)malloc(sizeof(double)*nd);

    //Startwerte initialisieren:
    x[0] = x_start;
    u[idx(0,nt,0, nd)] = 0.0;
    x[nd-1] = x_end;
    u[idx(0,nt, nd-1, nd)] = 0.0;
    for(int i=1; i<nd-1; i++){
        x[i] = x_start+i*h;
        
        u[idx(0,nt,i, nd)]= (double)(-N*(N+1))/pow(cosh(x[i]), 2);
    }
    t[0] = t_start;
    //printf("Initialisierung geschafft.\n");

    //ersten Zeitschritt mit anderen Differenzenquotienten durchführen.
    t[1] = t_start + d;
    for(int j=0; j<nd; j++){
        fak1 = 6.0*u[idx(0,nt, j,nd)]*(u[idx(0,nt,j+1, nd)]-u[idx(0,nt, j-1, nd)])/(2.0*h);
        fak2 = -(u[idx(0,nt, j+2, nd)] - 2.0*u[idx(0,nt, j+2,nd)] + 2.0*u[idx(0,nt, j-1,nd)] - u[idx(0,nt, j-2, nd)])/(2.0*h*h*h);
        u[idx(1,nt, j, nd)] = d*(fak1+fak2) + u[idx(0,nt, j, nd)];
    }
    //printf("Erster Schritt geschafft\n");
    

    for(int n=1; n<nt-1; n++){
        t[n] = t_start + n*d;
        for(int j=0; j<nd; j++){
            //Brechnung der neuen Punktes u_n^j
            //Setze den Rand immer auf 0
            if(j==0 || j==nd-1){
                u[idx(n+1,nt, j,nd)] = 0.0;
            }
            else{
                fak1 = ((u[idx(n,nt, j+1,nd)] + u[idx(n,nt,j,nd)] + u[idx(n,nt,j-1, nd)])/3.0) * ((u[idx(n,nt, j+1,nd)] - u[idx(n,nt,j-1, nd)])/(2.0*h));
                fak2 = (u[idx(n,nt,j+2,nd)] - 2.0*u[idx(n,nt,j+1,nd)] + 2.0*u[idx(n,nt,j-1,nd)] - u[idx(n,nt, j-2, nd)])/(2.0*h*h*h);
                u[idx(n+1,nt, j,nd)] = 2.0*d*(6.0*fak1 - fak2) + u[idx(n-1, nt,j, nd)];
                
            }
        }
    }

    for(int i=0; i<nd; i++){
        func_end[i] = u[idx(nt-1,nt,  i, nd)];
    }

    print_data_table_double(filename1, x, func_end, nd);
    if(filename2 !=NULL){
        print_table_1dim(filename2, u, nd, nt);
        print_data2file("ort.txt", x, nd);
        print_data2file("zeit.txt", t, nt);
    }
    //print_data2file(filename2, x, nd);
    //print_data2file(filename3, t, nt);
    free(u); free(x); free(t);free(func_end);
}

int main(){
    double x_start, x_end, t_start, t_end;
    double d,h;
    x_start = -30.0;
    x_end = 30.0;
    t_start = 0.0; 
    t_end = 1.0;
    printf("[main] Berechnung der numerischen Lösung für N=1 und N=2\n");
    printf("[main] Bitte wähle Parameter d und h mit d;h:");
    scanf("%lf;%lf", &d,&h);
    
    solve_kdv(t_start, t_end,d, x_start, x_end,h, 2, "funktion2.txt", NULL);
    solve_kdv(t_start, t_end,d, x_start, x_end,h, 1, "funktion1.txt", NULL);
    
    printf("[main] Berechnung der numerischen Lösung für N=3\n");
    printf("[main] Bitte wähle Parameter d und h mit d;h:");
    scanf("%lf;%lf", &d,&h);

    x_start = -5.0;
    x_end = 45.0;
    solve_kdv(t_start, t_end, d,  x_start, x_end, h,  3, "funktion3.txt", "verlauf_funktion3.txt");
    
    printf("[main] Berechnung der Funktionswerte für nicht-ganzzahlige N\n");
    x_start = -30.0;
    x_end = 30.0;
    t_start = 0.0; 
    t_end = 1.0;
    h=0.05;
    d=0.00003;
    char* funktions[]= {"funktion1_1.txt", "funktion1_2.txt", "funktion1_3.txt", "funktion1_4.txt", "funktion1_5.txt", "funktion1_6.txt", "funktion1_7.txt", "funktion1_8.txt", "funktion1_9.txt"};
    for(int i=1; i<10; i++){
        solve_kdv(t_start, t_end,d, x_start, x_end,h, 1.0+i*0.1, funktions[i-1], NULL);
    }
    
}