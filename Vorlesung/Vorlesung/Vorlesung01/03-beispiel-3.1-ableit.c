/* Datei: beispiel-3.1-ableit.c Datum: 16.4.2016 */  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>       /* mathematische Funktionen */ 

int main()   
{ double x=1.0,h=0.001,f0,f1,f2,dfa,dfb,dfc;      
  
  f0=sin(x-h);
  f1=sin(x);
  f2=sin(x+h);

  dfa=(f2-f1)/h;
  dfb=(f1-f0)/h;
  dfc=(f2-f0)/(2*h);
 
  /* gebe Differenzenquotienten (verschiedene Varianten aus), vergleiche mit exakter Ableitung 
     wähle Exponentialdarstellung für Ausgabe der Gleitkommazahlen  */

  printf("%15.6e   %15.6e   %15.6e   %15.6e \n",dfa,dfb,dfc,cos(x));
 
  return 0;
}

/* Ergebnis :   numerische Ableitung von sin(x) mit Methode (a),(b) und (c)

       5.398815e-01        5.407230e-01        5.403022e-01        5.403023e-01

*/
