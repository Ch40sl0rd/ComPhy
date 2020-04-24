/* Datei: beispiel-3.3-secant.c  Datum: 19.4.2020 */  
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Routine fuer Nullstellensuche mit dem Sekantenverfahren */ 
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
  
double f(double x)      /* eine Beispielfunktion */ 
{                       /* beide sind double Funktionen */ 
    return x*log(x)-x;  /* mit einem double Parameter */
}


int main()
{ 
  
  double a,b;              /* fuer die Startwerte */ 
  double exact,diff,res;  /* fuer die Ergebnisse */ 
  int n;                  /* Anzahl der Schritte */ 


  printf("Bitte geben Sie a,b ein: ");
  scanf("%lf %lf",&a,&b); 

  printf("       sekant         exakt          diff\n\n");

  /* Suche Nullstelle und speichere Ergebnis in "res" 
     Referenzdefinitionen bei "secant": 
       Es werden automatisch die Adressen der Objekte f und n uebergeben  
       n wird veraendert und enthaelt die Anzahl der Schritte nach Aufruf !!! */ 

  res=secant(a,b,&f,&n);       
                         
  exact=exp(1.0);          /* Vergleich mit exaktem Ergebnis */
  diff=fabs(res-exact);    
  
  printf("%15d   %15.6e   %15.6e   %15.6e \n",n,res,exact,diff);

  return 0;
}
/* Ergebnis: 
Bitte geben Sie a,b ein: 1 2 
         sekant          exakt           diff

    8 2.71828183e+00 2.71828183e+00 4.44089210e-16
*/
