/* Datei: beispiel-4.1-euler.c    Datum: 23.4.2012 */  

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<limits.h>

/* Routine, die rechte Seite der Dgl definiert. 
   neq: Anzahl der Gleichungen
   t : "Zeit" zur der rechte Seite benoetigt ist
   y : Loesung y zu dieser "Zeit" (Feld der Laenge neq)
   f : Ergebnis fuer die rechte Seite (Ausgabe) (Feld der Laenge neq)
*/ 

void derhs(int neq,double t,double *y,double *f)
{
  /* sicherstellen, dass Anzahl der Gleichungssystem wie erwartet ist */ 
  if(neq!=1) 
    {
      printf("neq passt nicht!\n");
      abort();
    } 

  /* es wird angenommen, dass die aufrufende Funktion den Speicherplatz f 
     bereitstellt 
     hier   dy/dt = - t * y(t) */ 

  f[0]=-t*y[0];

}

/* Schritt nach dem Euler-Verfahren 
   neq: Anzahl der Gleichungen
   h   : Schrittweite
   t   : "Zeit" beim Start 
   y   :  Loesung bei t ,  (Aufruf) 
          Loesung bei t+h, (Rueckgabe)
          Feld der Laenge neq  
   f   : Hilfsfeld der Laenge neq
   derhs : Zeiger auf Funktion, die rechte Seite bestimmt
*/   
   
void euler(int neq, double h, double t, double *y,double *f, 
           void (*derhs) (int, double ,double *, double*))
{
  int i;
  /* rechte Seite bestimmen f(t,y(t)) und bei f speichern */
 
  (*derhs)(neq,t,y,f);

  /* Schritt einfaches Eulerverfahren:
      y(t+h) = y(t) + h * f(t,y(t)) */ 

  for(i=0;i<neq;i++)
    {
      y[i]+=h*f[i];     
    }

}

int main()
{
  double h,t0,y0,tend,tstep;  /* Schrittweite, startpunkt, Startwert
                                 Endpunkt, Schritt fuer Ausgabe */ 
  int neq=1;                  /* feste Vorgabe der Anzahl der Gleichungen */ 
  double exact,diff;    /* Variablen, um Ergebnis zu speichern und vergleichen*/  
  double *y,*f;         /* Zeiger auf Speicherplaetze, die double enthalten */

  double t,tprint,eps=1.0E-4; 

  /* Eingabe der Parameter */ 
  printf("Bitte geben Sie h,t0,y0,tend und tstep ein: \n");
  scanf(" %le %le %le %le %le",&h,&t0,&y0,&tend,&tstep);  

  /* malloc reserviert Speicher fuer Feld mit y Werten und Hilfsfeld */ 
  y=malloc(sizeof(double)*neq);     
  f=malloc(sizeof(double)*neq);     

  printf("\n   %20s %20s %20s %20s \n","t","exact","dgl","diff");

  
  y[0]=y0;
  tprint=t0;

  for(t=t0;t<=tend;t+=h)
    {
      if(t-tprint>=-eps) /* Naechsten Ausgabepunkt erreicht?*/ 
	{
          exact=exp(-t*t*0.5);           /* known exact value */ 
          diff=fabs(exact-y[0])/exact;   /* and rel. error */
	  printf("   %20.5le %20.5le %20.5le %20.5le \n",t,exact,y[0],diff);
	  tprint+=tstep;  /* Ausgabe und naechsten Punkt bestimmen */
	}
      
      /* Dgl.schritt ausfuehren */       
      euler(neq,h,t,y,f,&(derhs));

    }    

  
  free(y);
  free(f);

}
