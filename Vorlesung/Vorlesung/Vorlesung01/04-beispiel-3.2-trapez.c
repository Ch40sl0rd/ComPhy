/* Datei: beispiel-3.2-trapez.c   Datum: 19.4.2020 */  
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Definition der zu integrierenden Funktion */

double f(double x)
{ 
  return exp(x);
}

/* Routine, die Gitterpunkte und Gewichte fuer ein Intervall [a,b] festlegt */

void trapez(int n, double a, double b, double *xp, double *wp)
/* n legt die Anzahl der Stuetzstellen fest. 
   a,b sind "normale" double Parameter
   xp und wp sind Zeiger(Pointer) auf ein double Feld, 
   es wird die Adresse des Feldes gespeichert !!!  */
{ 
  int i;
  double h;
  
  h=(b-a)/(double)(n-1);    /* Berechne Schrittweite */ 

  for(i=1;i<n-1;i++)
    { 
      xp[i]=a+i*h;      /* xp = Anfangsadresse des Feldes */     
      wp[i]=h;          /* xp[i] = Nehme die Speicherstelle, */
                        /*      die i Speicherstellen weiter liegt */	  	  
    }

  xp[0]=a;          /* Lege Punkte und Gewichte am Rand fest */
  wp[0]=h/2.0;

  xp[n-1]=b; 
  wp[n-1]=h/2.0; 

}

int main()
{ 
  double a,b;             /* Intervallgrenzen */
  int  i,n;               /* Schleifenvariable,Anzahl der Stuetzstellen */
  double exact,diff,sum;  /* Variablen, um Ergebnis zu speichern */  
  double *x,*w;         /* Zeiger auf Speicherplaetze, die double enthalten */

  printf("Bitte geben Sie a,b und n ein: ");  /* Eingabe der Parameter */
  scanf("%lf %lf  %d",&a,&b,&n);

  
  
  x=(double *)malloc(n*sizeof(double));     /* Anzahl der Stuetzstellen bestimmt Laenge des Feldes */
  w=(double *)malloc(n*sizeof(double));     /* waehrend des Programmlaufes ! */ 
                                            /* malloc erledigt diese Aufgabe ("allozieren") */

  trapez(n,a,b,x,w);   /* uebergibt n = Anzahl der Punkte a,b Intervallgrenzen */ 
                       /* Adressen der mit new allozierten Speicherbereiche */ 
  
  /* Gitterpunkte und Gewichte sind nun bei x und w gespeichert */ 
  /* Diese kann man fuer beliebige Funktionen benutzen */ 

  sum=0.0;              /* Bestimmung eines Integrals mit den Gitter und Gewichten */

  for(i=0;i<n;i++)
   {
    sum+=f(x[i])*w[i];   /* "+=" Operator summiert f(xi)*w(xi) auf Summe auf */
   }

  exact=exp(b)-exp(a);
  diff=fabs(sum-exact);

  printf("N      trapez      exact       diff \n\n");
  printf("%d      %15.6e      %15.6e       %15.6e \n",n,sum,exact,diff);

  
  free(x);   /* Speicherbereich wieder freigeben */ 
  free(w);   /*  ("deallozieren") */

  return 0;
   
}

/* Ergebnis:  fuer a b n = 0.0 1.0 30
                   N              trapez               exact                diff

                  30    1.7184520869e+00    1.7182818285e+00    1.7025840042e-04
*/
