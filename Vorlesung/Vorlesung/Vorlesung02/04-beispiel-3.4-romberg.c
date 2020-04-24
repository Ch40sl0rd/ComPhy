/* Datei: beispiel-3.4-romberg.c    Datum: 19.4.2020 */  
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int mmax=10;  /* definiere globale Variable mit maximaler Anzahl der Schritte */

/* Definition der zu integrierenden Funktion */

double f(double x)
{ 
  return exp(x);
}

/* Funktion die Trapezsumme mit n Funktionspunkten bestimmt fuer Integralgrenzen a,b und Funktion func */

double T(int n,double a, double b, double (*func)(double))
{ double sum,h;
  int i;

  h=(b-a)/(n-1);

  sum=0.5*func(a)+0.5*func(b);

  for(i=1;i<n-1;i++)
    { sum+=func(a+h*i); }
      
  return h*sum; }  

/* Funktion, die eindeutigen Index definiert von ik -> index_ik  fuer  i>=k und i,k<=mmax */ 
int index_ik(int i,int k)
   { return (i-k + k*(mmax+1)-k*(k-1)/2); }


/* Routine, die Romberg Integration bis zu einer 
   "Genauigkeit" eps durchführt. 
   n0 ist bei Aufruf die  Anzahl der Stuetzstellen im ersten Schritt 
   und beim Verlassen die maximale Anzahl der benutzten Stuetzstellen
   a, b die Integralgrenzen und f die zu integrierend Funktion    */

double romberg(int *n0, double a, double b, double (*func)(double),double eps)
{ int k,m,n;  /* fuer Indizes */ 
  double *h;  /* Schrittweiten fuer j=0,...,mmax */  
  double *Tsum; /* Trapezsummen fuer j=0,...,mmax */  
  double *tildeT; /* Neville-Schema tilde T_{jk}  bei h=0 */
  double result;
 
  h=(double *)malloc((mmax+1)*sizeof(double));       /* Speicher fuer h in Schritt m */
  Tsum=(double *)malloc((mmax+1)*sizeof(double));    /* und T fuer diese h */
  tildeT=(double *)malloc(((mmax+1)*(mmax+2))/2*sizeof(double));  /* Speicher fuer Neville Schema */ 

  h[0]=(b-a)/(double)(*n0-1);
  n=*n0;
  Tsum[0]=T(n,a,b,func);
  tildeT[index_ik(0,0)]=Tsum[0];

  for(m=1;m<=mmax;m++) 
    { h[m]=h[m-1]/2;         /* Trapezsumme fuer halbiertes h */ 
      n=2*n;
      Tsum[m]=T(n,a,b,func);

      tildeT[index_ik(m,0)]=Tsum[m]; /* fuer  i=m und k=0 */ 

      /* generate tildeT i=m k=1,...,m */ 
      for(k=1;k<=m;k++)
	{ tildeT[index_ik(m,k)]=-h[m-k]*h[m-k]/(h[m]*h[m]-h[m-k]*h[m-k])*tildeT[index_ik(m,k-1)]
	                  + h[m]  *h[m]  /(h[m]*h[m]-h[m-k]*h[m-k])*tildeT[index_ik(m-1,k-1)];
         }

      printf("%5d   %15.6e  %15.6e \n",m,tildeT[index_ik(m,m)],Tsum[m]);  /* just to observe convergence */
      if(fabs(tildeT[index_ik(m,m)]-tildeT[index_ik(m-1,m-1)])<=eps) break; /* stop when accuracy reached */

    }
  *n0=n;

  result=tildeT[index_ik(m,m)];

  free(tildeT);
  free(Tsum);
  free(h);

  return result;

}

int main()
{ 
  double a,b;             /* Intervallgrenzen */
  int n;                  /* Stuetzstellen am Anfang */
  double exact,diff,sum;  /* Variablen, um Ergebnis zu speichern */  

  printf("Bitte geben Sie a,b und n ein: ");  /* Eingabe der Parameter */
  scanf("%lf %lf  %d",&a,&b,&n);


  sum=romberg(&n,a,b,&f,1.0E-6);   /* uebergibt n = Anzahl der Punkte beim Start,  a,b Intervallgrenzen */ 
                                /* Adressen der Funktion und Genauigkeit */ 
  
  exact=exp(b)-exp(a);
  diff=fabs(sum-exact);

  printf("Nmax      romberg       exact       diff \n\n");
  printf("%d        %15.6e      %15.6e       %15.6e \n",n,sum,exact,diff);

  return 0;   
}
/* Ergebnis: 
 Bitte geben Sie a,b und n ein: 0 1 2 
    1      1.692503e+00     1.734162e+00 
    2      1.718509e+00     1.721203e+00 
    3      1.718237e+00     1.718918e+00 
    4      1.718277e+00     1.718431e+00 
    5      1.718281e+00     1.718318e+00 
    6      1.718282e+00     1.718291e+00 
Nmax      romberg       exact       diff 

128           1.718282e+00         1.718282e+00          8.316105e-08 (Trapezregel benoetigt etwa 1200 Punkte)*/
