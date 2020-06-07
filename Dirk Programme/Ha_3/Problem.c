

//Numerov Verfahren zu Aufgabe 8

//Einfuegen der Standardbibliotheken und Pakete
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Define von M_PI wegen der nicht vorhandenen Math_PI in Visual studios
#ifndef M_PI
#define M_PI 3.141592653589793238
#endif

// C-Programme haben Probleme in MS VISUAL deswegen hier die Fehlercodeblockierung
#pragma warning(disable : 4996)

//Im folgenden Programm wurde die numerov Methode aus der Vorlesung 8 entnommen

/* beispiel-4.4-numerov.c   Datum:  06.05.2020 */
/* dieser Code implementiert ein vereinfachtes Numerov Problem,
   Erweiterung notwendig bei "Korrektur notwendig fuer g!=0 !!!!" */

//####################################################################################################################


double secant(double x1, double x2, double (*func)(double), int* schritt)
/* x1,x2     Startwerte
   func      ist die "Referenz" auf eine Funktion mit einem double Parameter,
             die double zurueckgibt (Referenz = Adresse der Funktion)
   schritt   ist auch Referenz: Veraenderungen an der Variable wirken sich auf
                                das aufrufende Programm aus !!! */
{
    const double tol = 1e-12; /* geforderte Genauigkeit als Konstante */
    double xn;              /* neuer Schaetzwert */

    *schritt = 0;  /* noch kein Schritt */

    do
    {         /* naechster Schaetzwert x1,x2 -> xn*/
        xn = x2 - func(x2) * (x2 - x1) / (func(x2) - func(x1));
        x1 = x2;     /* bereite den naechsten Schritt vor:  x2 -> x1 */
        x2 = xn;     /*                                     xn -> x2 */

        (*schritt)++;      /* Schritte=Schritte+1 */
    } while (fabs(x2 - x1) > tol);   /* solange Genauigkeitsziel nicht erreicht */

    return xn;   /* Gebe Nullstelle zurueck */
}


/* Code fuer die Loesung der Poisson-Gleichung mit gaussfoermiger Ladungsverteilung */

/* wegen Nullstellensuche definiere globale Felder, die in
   Funktion zur Nullstellensuche genutzt werden koennen */

//double r0bound;                               /* definiert radius der homogenen Ladungsverteilung */
double* r_array, * g_array, * s_array, * y_loesung; /* Zeiger auf Felder fuer Numerov und Loesung */
double schrittweite_h;                        /* benutzte Schrittweite */
int num_r;                                    /* Anzahl der Stuetzstellen */

//Funktionen f�r numerov Verfahren der Aufgabe 8 #####################################################################
//Zu loesende Gleichung die dem V(x) aus der Aufgabenbeschreibung entspricht
double gfunc(double x)
/* zu loesende Gleichung: y''(r)+g(r)*y(r)=s(r)
   hier Definition der Funktion g(r)   */

{
    return 60 * pow(cos(M_PI * x), 16);
}

//8.3 zu l�sende Gleichung mit �u�erem Elektrischen Feld
double gfuncel(double x, double epsilon)
/* zu loesende Gleichung: y''(r)+g(r)*y(r)=s(r)
   hier Definition der Funktion g(r)   */

{
    return 60 * pow(cos(M_PI * x), 16)+ epsilon*x;
}
//####################################################################################################################

double sfunc(double r)
// zu loesende Gleichung: y''(r)+g(r)*y(r)=s(r)
// hier Definition der Funktion s(r)
{

    return 0;
}


/* fuer konstante Verteilung, vorsicht Unstetigkeit fuehrt zu Ungenauigkeit !!!!
double sfunc(double r)
// zu loesende Gleichung: y''(r)+g(r)*y(r)=s(r)
//   hier Definition der Funktion s(r)
// Beispiel einer unstetigen Funktion
  {
    if(r>r0bound)
      {
    return 0.0;
      }
    else
      {
    return -4.0*r*1.0/(4.0*M_PI/3.0*r0bound*r0bound*r0bound);
      }
      } */

double init_numerov(double a, double b, int n, double* r, double* g, double* s)
/* Funktionen bereitet Anwendung des Numerov Verfahrens vor
   a,b,n sind Intervallgrenzen und Anzahl der Stuetzstellen
   bei Verlassen der Routine beinhalten die Felder r,g und s
   Stuetzstellen, Funktion g an den Stuetzstellen und s an den Stuetzstellen
   Rueckgabewert ist die Schrittweite h */
{
    int i;
    double h;

    h = (b - a) / (n - 1);                /* bestimme Intervalllaenge */
    for (i = 0; i < n; i++)            /* belege die Felder mit den Werten fuer r,g und s*/
    {
        r[i] = a + i * h;
        g[i] = gfunc(r[i]);
        s[i] = sfunc(r[i]);
    }

    return h;  /* Rueckgabewert ist Schrittweite h */
}

double init_numerovel(double a, double b, int n, double* r, double* g, double* s, double epsilon)
/* Funktionen bereitet Anwendung des Numerov Verfahrens vor
   a,b,n sind Intervallgrenzen und Anzahl der Stuetzstellen
   bei Verlassen der Routine beinhalten die Felder r,g und s
   Stuetzstellen, Funktion g an den Stuetzstellen und s an den Stuetzstellen
   Rueckgabewert ist die Schrittweite h */
    // Hier wird das externe Elektrische Feld �usgedr�ckt durch epsilon beachtet
{
    int i;
    double h;

    h = (b - a) / (n - 1);                /* bestimme Intervalllaenge */
    for (i = 0; i < n; i++)            /* belege die Felder mit den Werten fuer r,g und s*/
    {
        r[i] = a + i * h;
        g[i] = gfuncel(r[i], epsilon);
        s[i] = sfunc(r[i]);
    }

    return h;  /* Rueckgabewert ist Schrittweite h */
}


void numerovup(double* r, double* g, double* s, double h, int n, int steps, double y0, double y1, double* y)
/* Funktion nutzt das Numerov Verfahren um fuer gegebene Stuetzstellen r
   und Funktionen g und s die Loesung y zu finden.
   r, g und s sollten mit init_numerov vorbereitet werden
   h ist die Schrittweite (auch aus init_numerov)
   n ist die Anzahl der Stuetzstellen
   steps legt fest, wieviele Numerovschritte ausgefuehrt werden
   y0 und y1 sind die Startwerte y[0] und y[1], die Routine legt dann y[2] .. y[steps+1] fest
   y ist ein Feld der Laenge  steps+2, das die Loesung bei Rueckkehr enthaelt
*/
{
    int i;
    double fakt_u_np1, fakt_u_n, fakt_u_nm1; /* Variablen fuer Numerov Faktoren */
    double fakt_s;

    y[0] = y0;  /* belege erste Funktionswerte mit den Startwerten  */
    y[1] = y1;

    for (i = 1; i < steps + 1; i++)   /* betrachte in den Schritten y(i-1) und y(i) um y(i+1) zu berechnen */
    {
        fakt_u_np1 = 1.0;     /* Faktor bei y(i+1),  Korrektur notwendig fuer g!=0 !!!! */
        fakt_u_n = 1.0;       /* Faktor bei y(i),    Korrektur notwendig fuer g!=0 !!!! */
        fakt_u_nm1 = 1.0;     /* Faktor bei y(i-1)   Korrektur notwendig fuer g!=0 !!!! */
        fakt_s = h * h / 12.0 * (s[i + 1] + 10.0 * s[i] + s[i - 1]);   /* Faktor mit s */

        y[i + 1] = (fakt_s + 2.0 * fakt_u_n * y[i] - fakt_u_nm1 * y[i - 1]) / fakt_u_np1; /* Vorwaertsiteration */
    }

}

void numerovdown(double* r, double* g, double* s, double h, int n, int steps, double yn, double ynm1, double* y, double lambda)
/* Funktion nutzt das Numerov Verfahren um fuer gegebene Stuetzstellen r
   und Funktionen g und s die Loesung y zu finden.
   r, g und s sollten mit init_numerov vorbereitet werden
   h ist die Schrittweite (auch aus init_numerov)
   n ist die Anzahl der Stuetzstellen
   steps legt fest, wieviele Numerovschritte ausgefuehrt werden
   yn und ynm1 sind die Startwerte y[n-1] und y[n-2], die Routine legt dann y[n-3] .. y[n-steps-2] fest
   y ist ein Feld der L�nge n, das die Loesungen
*/
{
    int i;
    double fakt_u_np1, fakt_u_n, fakt_u_nm1; /* Variablen fuer Numerov Faktoren */
    double fakt_s;

    y[n - 1] = yn;  /* belege erste Funktionswerte mit den Startwerten  */
    y[n - 2] = ynm1;

    for (i = n - 2; i > n - steps - 2; i--)   /* betrachte in den Schritten y(i-1) und y(i) um y(i+1) zu berechnen */
    {
        fakt_u_np1 = (1.0 + h * h / 12 * (2 * (-g[i + 1] + lambda)));     /* Faktor bei y(i+1),  Korrektur notwendig fuer g!=0 !!!! */
        fakt_u_n = (1.0 - 5 * h * h / 12 * (2 * (-g[i] + lambda)));       /* Faktor bei y(i),    Korrektur notwendig fuer g!=0 !!!! */
        fakt_u_nm1 = (1.0 + h * h / 12 * (2 * (-g[i - 1] + lambda)));     /* Faktor bei y(i-1)   Korrektur notwendig fuer g!=0 !!!! */
        fakt_s = h * h / 12.0 * (s[i + 1] + 10.0 * s[i] + s[i - 1]);   /* Faktor mit s */

        y[i - 1] = (fakt_s + 2.0 * fakt_u_n * y[i] - fakt_u_np1 * y[i + 1]) / fakt_u_nm1; /* Rueckwaertsiteration */
    }

}

double f_randbed(double para,double lambda)
/* Nullstelle dieser Funktion signalisiert, dass die Loesung
   konsistent mit der Randbedinung ist
   para ist der freie Parameter der Loesung
   hier: para ist der Konstantewert der Loesung bei gro�en r
   Die Funktion nutzt globale Variabeln, um die Felder
   r, g, s, h und n zu erhalten
   und feld y_loesung fuer die Loesung */
{
    numerovdown(r_array, g_array, s_array, schrittweite_h, num_r, num_r - 2, 0, para, y_loesung, lambda);
    return y_loesung[0];  /* zweite Randbedingung ist y[0] = 0 */
}
//Funktion nach b umgestellt und 1/ wurzel(Summe der Integrale) gerechnet
double Normfakt(double a, double b, double* y) {
    double c = (b - a) / num_r;
    double sum = 0;
    for (int i = 0; i < num_r; i++){
        sum = sum + c * (pow(y[i], 2));
    }

    return 1 / sqrt(sum);

}

int main()
{
    double rmax, para, zero;  /* maximales r und konsistente Randbedingung und Null aus Auswertung der Randbedingungsfunktion */
    int i, num_schritt; /* Anzahl der Sekantenverfahren Schritte */
    double lambda = 0;
    int ort = 0;

    //double rho;  /* fuer Ausgabe der Ladungsverteilung */

    //printf("# Geben Sie rmax,r0 und die Anzahl der Stuetzstellen ein:");
    //scanf(" %le %le %d", &rmax, &r0bound, &num_r);


    printf("Zu welchem Aufgabenteil m�chten sie 8.(1,2,3)?\n");
    printf("Oder moechten Sie beenden ? (4)\n");
    scanf("%d", &ort);

    r_array = malloc(num_r * sizeof(double));   /* alloziere Speicher fuer r, g und s */
    g_array = malloc(num_r * sizeof(double));
    s_array = malloc(num_r * sizeof(double));
    y_loesung = malloc(num_r * sizeof(double));

    //schrittweite_h = init_numerov(0.0, rmax, num_r, r_array, g_array, s_array);  /* belege r,g und s mit Werten*/

    //Aufgabe 8.1 ####################################################################################################
    if (ort == 1)
{
    //Länge=8
    rmax = 8;
    schrittweite_h = init_numerov(0.0, rmax, num_r, r_array, g_array, s_array);
    int o = 1;
    //Bestimmen der Energieeigenwerte im Bereich 0 bis maxV_per(x) = 60 in Schrittweite von 0.01

    for (double i = 0; i < 60; i = i + 0.01) {

    f_randbed(1 * 10e-6, i);

        if (fabs(y_loesung[0]) <= 10e-7) {

        printf("Energieeigenwert %d: %lf\n", o, i);

        printf("Die Normierungskonstante beträgt: N=%10e\n", Normfakt(0,8,y_loesung));

        o++;
        }

    }



}

//Aufgabe 8.2 ####################################################################################################
else if (ort == 2)
{
    //Längeneingabe 16,32....
    printf("Bitte neue Länge eingeben(16,32...):");
    scanf("%lf", &rmax);
    schrittweite_h = init_numerov(0.0, rmax, num_r, r_array, g_array, s_array);
    int o = 1;
    for (double i = 0; i < 60; i = i + 0.01) {

        f_randbed(10e-6, i);

        if (fabs(y_loesung[0]) <= 10e-7) {

        printf("Energieeigenwert %d: %lf\n", o, i);

        o++;
        }

    }

}



//Aufgabe 8.3 ####################################################################################################
else if (ort == 3)
{
    rmax = 8;
    int o = 1;
    double epsilon = 0;
    int Band = 0;
    double zwew = 0;
    while (epsilon < 0.5){
        epsilon = epsilon + 0.001;
        schrittweite_h = init_numerovel(0.0, rmax, num_r, r_array, g_array, s_array, epsilon);

        for (double i = 0; i < 60 && Band < 2; i = i + 0.01) {

            f_randbed(10e-6, i);

            if (fabs(y_loesung[0]) <= 10e-7) {

                if (Band<2 && zwew != 0 && fabs(i - zwew))
                {
                    printf("Energieeigenwert %lf: %lf,\n", zwew, i, fabs(i - zwew));

                    Band++;
                }

                zwew = i;
                o++;

            }

        }

    }

}


else
{
    printf("Sie haben sich verschrieben");

}

    free(r_array);
    free(g_array);
    free(s_array);
    free(y_loesung);

    return 0;
}
