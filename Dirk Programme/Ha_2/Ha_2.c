// gcc Ha_2.c libh2_cip.a -o Ha_2 -lm
// ./Ha_2
// Autor: Dirk Knott


//Einfuegen der Standardbibliotheken und Pakete

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<float.h>
#include<limits.h>
#include<complex.h>
#include "h2.h"
// C-Programme haben Probleme in MS VISUAL deswegen hier die Fehlercodeblockierung
#pragma warning(disable : 4996)

//Entnommene Funktionen aus:
/* Datei: beispiel-4.2-rungekutta-2.c  Datum: 24.4.2012
   04-Runge-Kutta-4.c
*/

//####################################################################################################################

//H.4 Gedämpfte Oszillation

//H.4.2: Implementation der Simulation des Systems mit Hilfe des RK4-Verfahrens

/* Routine, die rechte Seite der Dgl definiert.
   neq: Anzahl der Gleichungen
   t : "Zeit" zur der rechte Seite benoetigt ist
   y : Loesung y zu dieser "Zeit" (Feld der Laenge neq)
   f : Ergebnis fuer die rechte Seite (Ausgabe) (Feld der Laenge neq)
*/
//Funktion mit bestimmten Gleichungssystem der zu lösenden DGL
void derhs(int neq, double t, double* y, double* f)
{
    /* sicherstellen, dass Anzahl der Gleichungssystem wie erwartet ist */
    if (neq != 2)
    {
        printf("neq passt nicht!\n");
        abort();
    }

    /* es wird angenommen, dass die aufrufende Funktion den Speicherplatz f
       bereitstellt
       Definition des bestimmten Glgsys. der DGL 1.Ordnung eines gedämpften harmon Oszi
       wie in Aufgabe vorgegeben d^2x/dt^2 = -w_0^2*x-2*g*dx/dt
       w_0=1 und gamma = 0.1 -> 2*gamma=0.2
       */

    f[0] = y[1];
    f[1] = -y[0]-0.2*y[1];
}



//#####################################################################################################################

//H.5:Kopplung an ein externes Feld
//H.5.1:Implementation der Simulation des Systems für den Spezialfall
//f(t)= alpha * cos(omega * t)

void derhs_2(int neq, double t, double* y, double* f)
{
    /* sicherstellen, dass Anzahl der Gleichungssystem wie erwartet ist */
    if (neq != 2)
    {
        printf("neq passt nicht!\n");
        abort();
    }

    /* es wird angenommen, dass die aufrufende Funktion den Speicherplatz f
       bereitstellt
       Definition des bestimmten Glgsys. der DGL 1.Ordnung eines gedämpften harmon Oszi
       wie in Aufgabe vorgegeben d^2x/dt^2 = -w_0^2*x-2*g*dx/dt
       w_0=1 und gamma = 0.1 -> 2*gamma=0.2
       */

    f[0] = y[1];
    f[1] = -y[0] - 0.2 * y[1]+0.5*cos(1.5*t);
}

//########################################################################################################################
//H.6:Tune in
//H.6.1:Erweitern der Simulation
//Erweitern der Simulation mit fextern eingefügt
void exter(int neq, double t, double* y, double* f){
    /* sicherstellen, dass Anzahl der Gleichungssystem wie erwartet ist */
    if (neq != 2)
    {
        printf("neq passt nicht!\n");
        abort();
    }

    /* es wird angenommen, dass die aufrufende Funktion den Speicherplatz f
       bereitstellt
       Definition des bestimmten Glgsys. der DGL 1.Ordnung eines gedämpften harmon Oszi
       wie in Aufgabe vorgegeben d^2x/dt^2 = -w_0^2*x-2*g*dx/dt
       w_0=1 und gamma = 0.1 -> 2*gamma=0.2
       */

    f[0] = y[1];
    f[1] = -y[0] -0.2 * f[0] + fext(t);
}



//########################################################################################################################

/*Runge-Kutta-Verfahren 4. Ordnung*/
void RuKu_4(int nDifEqu, /* # der Differentialgleichungen */
    double h,    /* Schrittweite */
    double t,    /* Kurvenparameter */
    double y[],  /* Bahnkurve [nDifEqu] */
                 /* Eingabe: y(t) */
                 /* Rueckgabe; y(t+h) */
    double yh[], /* Hilfsfeld [nDifEqu] */
    /* Hilfsfelder [nDifEqu]: */
    double k1[], double k2[], double k3[], double k4[],
    /* Funktion zur Berechnung der rechten Seite: */
    void (*derhs) (int, double, double[], double[])
) {

    /* Berechnet 1 Runge-Kutta Schritt (4. Ordnung) zur Loesung des DG-Systems:
       y'(t) = f(y(t),t); y \in R^n */

       /* Variablen */
    double h2;
    int i;

    h2 = 0.5 * h;

    //Berechnung der K in verschiedenen Ordnungen exakte Formeln siehe Tex
    (*derhs)(nDifEqu, t, y, k1);
    for (i = 0; i < nDifEqu; i++) {
        yh[i] = y[i] + h2 * k1[i];
    }
    (*derhs)(nDifEqu, t + h2, yh, k2);
    for (i = 0; i < nDifEqu; i++) {
        yh[i] = y[i] + h2 * k2[i];
    }
    (*derhs)(nDifEqu, t + h2, yh, k3);
    for (i = 0; i < nDifEqu; i++) {
        yh[i] = y[i] + h * k3[i];
    }
    (*derhs)(nDifEqu, t + h, yh, k4);
    for (i = 0; i < nDifEqu; i++) {
        y[i] += (h2 * (k1[i] + k4[i]) + h * (k2[i] + k3[i])) / 3;
    }

}



//Ausgabe:
int main()
{
    //Deklaration der Variablen
    double h, t0, y0, tend, tstep;  /* Schrittweite, startpunkt, Startwert,Endpunkt, Schritt fuer Ausgabe */
    int neq = 2;                  /* feste Vorgabe der Anzahl der Gleichungen */
    double exact, diff;    /* Variablen, um Ergebnis zu speichern und vergleichen*/
    double* y, * f, * k1, * k2, * k3, * k4;         /* Zeiger auf Speicherplaetze, die double enthalten */
    double t, tprint, eps = 1.0E-4;
    int ort;
    
    /* malloc reserviert Speicher fuer Feld mit y Werten und Hilfsfeld */
    y = malloc(sizeof(double) * neq);
    k1 = malloc(sizeof(double) * neq);
    k2 = malloc(sizeof(double) * neq);
    k3 = malloc(sizeof(double) * neq);
    k4 = malloc(sizeof(double) * neq);
    f = malloc(sizeof(double) * neq);
    
    //
    FILE *data= fopen("Text.txt","w");
    //Startpunkt für goto funktion Für das Benutzer*innenfreundliche Interface
    start:
    
    //Fragt den/die Benutzer*in danach in welchen Aufgabenteil er möchte und
    printf("Bitte geben Sie an zu welchem Aufgabenteil sie moechten(4.2(1),5.1(2),6.1(3) oder moechten sie Verlassen(4)?\n");
    scanf("%d", &ort);

    //Anfangsbedingungen aus Aufgabenstellung
    //Startzeit
    t0 = 0;
    //Startwert x'
    y[0] = 1.5;
    //Startwert x'
    y[1] = 2.75;
    //Print Zeit
    tprint = t0;
    //Omega Teil 6
    int omg = 1;
    //
    h = 0.01;

    if (ort == 1) {
        //Beschreibung der Aufgabe 4
        printf("Sie befinden sich in Aufgabenteil 1: Gedaempfte Oszillation\n");
        printf("Hier simulieren wir das System mit Hilfe des Runge Kutta Verfahrens\n");

        /* Eingabe der Parameter */
        printf("Bitte geben Sie h, tend und tstep ein: \n");
        scanf(" %le %le %le", &h, &tend, &tstep);

        //Tabellenbeschriftung
        printf("Die Ergebnisse lauten\n");
        printf("\n   %20s %20s %20s %20s \n", "t", "exact", "dgl", "diff");
        

        //For-Schleife die die DGL berechnet und den Vergleich mit dem exakten Ergebnis ausführt
        for (t = t0;t <= tend;t += h)
        {
            exact = 1.5 * pow(2.71828182846, -0.1 * t) * cos(pow(1 - 0.01, 0.5) * t) + 2.914609 * pow(2.71828182846, -0.1 * t) * sin(pow(1 - 0.01, 0.5) * t);/* known exact value */
                     
            diff = fabs(exact - y[0]) / exact;   /* and rel. error */

            if (t - tprint >= -eps) /* Naechsten Ausgabepunkt erreicht?*/
            {
                printf("   %20.5le %20.5le %20.5le %20.5le \n", t, exact, y[0], diff);
                tprint += tstep;  /* Ausgabe und naechsten Punkt bestimmen */
            }
            //fprintf(data, "%lf %lf\n", t, diff);
            /* Dgl.schritt ausfuehren */
            RuKu_4(neq, h, t, y, f, k1, k2, k3, k4, &(derhs));

        }
        goto start;
    }
    
    
    else if (ort == 2) {
        //Beschreibung der Aufgabe 5
        printf("Sie befinden sich in Aufgabenteil 2: Kopplung an ein externes Feld\n");
        printf("Simulation eines Systems mit der Einwirkung eines externen Feldes \n");

        //Eingabe der Parameter und Tabellenbeschriftung
        printf("Bitte geben Sie tend und tstep ein: \n");
        scanf("%le %le", &tend, &tstep);
        printf("\n   %20s %20s %20s %20s \n", "t", "dgl", "exact", "diff");

        //For-Schleife die die DGL berechnet und den Vergleich mit dem exakten Ergebnis ausführt
        for (t = t0; t <= tend; t += h)
        {

            /* known exact value */
            exact = -0.378214826021 * cos(1.5 * t) + 0.090771558245 * sin(1.5 * t) + 1.87821482602 * exp(-0.1 * t) * cos(pow(0.99, 0.5) * t) + 2.81577841162 * exp(-0.1 * t) * sin(pow(0.99, 0.5) * t);


            /* and rel. error */
            diff = fabs(exact - y[0]) / exact;


            if (t - tprint >= -eps) /* Naechsten Ausgabepunkt erreicht?*/
            {

                printf("   %20lf %20lf %20lf       %15.8e\n", t, y[0], exact, diff);
                tprint += tstep;  /* Ausgabe und naechsten Punkt bestimmen */
            }
            //fprintf(data, "%lf %lf\n", t, diff);
            /* Dgl.schritt ausfuehren */
            RuKu_4(neq, h, t, y, f, k1, k2, k3, k4, &derhs_2);

        }


        goto start;
    }
    
    
    else if (ort == 3) {

        //Beschreibung der Aufgabe 6
        //Bemerkung 6.3 nicht enthalten
        printf("Sie befinden sich in Aufgabenteil 3: Tune In\n");
        printf("Simulation eines Systems mit der Einwirkung einer externen Kraft \n");
        

        //Eingabe der Parameter und Tabellenbeschriftung
        printf("Bitte geben Sie tend und tstep ein: \n");
        scanf("%le %le", &tend, &tstep);
        printf("\n   %20s %20s %20s %20s \n", "t", "dgl","x'","Leistung");
        double P = 0;
        h = 0.001;
        //Zählvariablen
        int i = 0;
        int ns = 1;
        int g = 0;
        //Allokation des Speichers für die Leistung und 
        double *psammeln = malloc(sizeof(double) * 100000);
        double *pstriche = malloc(sizeof(double) * 60);
        double vorherig = 0;








        //For-Schleife die die DGL berechnet und die Leistung berechnet
        for (t = t0; t <= tend; t += h)
        {   
            vorherig = y[0];
            P = y[1] * fext(t);
            psammeln[g] = P;
            g++;
            
            if (t - tprint >= -eps) /* Naechsten Ausgabepunkt erreicht?*/
            {

                printf("   %20lf %20lf %20lf %20lf\n", t, y[0],y[1],P);
                tprint += tstep;  /* Ausgabe und naechsten Punkt bestimmen */
            }
            //fprintf(data, "%lf %lf\n", t, diff);
            /* Dgl.schritt ausfuehren */
            RuKu_4(neq, h, t, y, f, k1, k2, k3, k4, &exter);
            //Abfrage nach Nulldurchgang
            if ((y[0] > 0 && vorherig < 0 && i == 1) || (y[0] < 0 && vorherig > 0 && i == 1))
            {
                ns++;//Hochzählen NST
            }
                // Abfrage Initialnst
            if ((y[0] > 0 && vorherig < 0 && i == 0) || (y[0] < 0 && vorherig > 0 && i == 0))
            {
                i = 1;
                g = 0;
            }
            if (ns == 3) {

                double temp = 0;
                for (int l = 0;l <= g; l++) {
                    temp = temp + psammeln[l];
                }
                printf("Mittelwert der Aktuellen Periode %lf \n", temp / g);   //Mittelwert berechnen
                //Zurücksetzen der Var
                g = 0;
                ns = 1;
                
            }
        }
        goto start;
    }


    else if (ort == 4) {
        //Aufruf in Interface zum Verlassen des Programms
        goto stopp;
    }
    else {
        printf("Sie haben sich vertippt!");

        goto start;
    }
    



stopp:
//fclose(data);
    // Freigabe des allokierten Speichers
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(y);
    free(f);
    //Verlassen
    return 0;


}