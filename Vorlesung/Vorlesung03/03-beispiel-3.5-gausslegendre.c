/* Datei: beispiel-3.5-gausslegendre.c   Datum: 24.4.2020 */  
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* GNU Scientific Library (gsl) Routinen erlauben Gauss-Legendre Punkte zu bestimmen */ 
#include <gsl/gsl_integration.h>

void gausslegendre(double a,double b,double *x,double *w,size_t n)
 { gsl_integration_glfixed_table *xwtable;
   size_t i;
   double xi,wi;

   xwtable=gsl_integration_glfixed_table_alloc(n);

   if(xwtable==NULL)
     {
      printf("Problem with Gauss-Legendre\n");
      abort();
     }

   for(i=0;i<n;i++)
    {
      gsl_integration_glfixed_point (a, b, i, &x[i], &w[i], xwtable);
    }
   
   gsl_integration_glfixed_table_free(xwtable);  
    
  }

/* Definition der zu integrierenden Funktion */

/* Definiere globale Variabel, um Parameter der Funktion festzulegen 
   diese Variablen sind von allen Funktionen zugaenglich!!! */

double aconst=0.5,bconst=1;  /* Werte sind bei Programmstart vorgegeben */


double f(double x)
{ 
  return aconst/(bconst+x*x);   
       /* ... auch hier kann ich globale Variablen benutzen */
}

int main()
{ 
  double a,b;                 /* Intervallgrenzen */
  int  i,n;                   /* Schleifenvariable,Anzahl der Stuetzstellen */
  double  exact,diff,sum;     /* Variablen, um Ergebnis zu speichern */  
  double *x,*w;               /* Zeiger auf Speicherplaetze, die double enthalten */

  printf("Bitte geben Sie a,b und n ein: ");  /* Eingabe der Parameter */ 
  scanf("%lf %lf %d",&a,&b,&n);

  x=(double *) malloc(n*sizeof(double));     /* Anzahl der Stuetzstellen bestimmt Laenge des Feldes */
  w=(double *) malloc(n*sizeof(double));     /* waehrend des Programmlaufes ! */ 

  gausslegendre(a,b,x,w,n);   /* uebergibt n = Anzahl der Punkte a,b Intervallgrenzen */ 
                       /* Adressen der mit new allozierten Speicherbereiche */ 
  
  /* Gitterpunkte und Gewichte sind nun bei x und w gespeichert */ 
  /* Diese kann man fuer beliebige Funktionen benutzen */ 

  /* Lege jetzt die Parameter der Funktion durch 
     Definition der globalen Variabeln fest */

  aconst=1.0;
  bconst=1.0;

 
  sum=0.0;              /* Bestimmung eines Integrals mit den Gitter und Gewichten */

  for(i=0;i<n;i++)
   {
    sum+=f(x[i])*w[i];   /* "+=" Operator summiert f(xi)*w(xi) auf Summe auf */
   }

  exact=atan(b)-atan(a);
  diff=fabs(sum-exact);
  

  printf("N      gauleg      exact       diff \n\n");
  printf("%d      %15.6e      %15.6e       %15.6e \n",n,sum,exact,diff);


  free(x);   /* Speicherbereich wieder freigeben */
  free(w);   /*  ("deallozieren") */
  
  return 0; 
}

/* Compilation auf cip Rechnern:    gcc beispiel-3.5.c -lgsl -lblas -o beispiel-3.5 
Ergebnis:  fuer a b n = 0.0 1.0 5
Bitte geben Sie a,b und n ein: 0 1 5
N      gauleg      exact       diff 

5         7.853982e-01         7.853982e-01          3.426266e-09 
*/


