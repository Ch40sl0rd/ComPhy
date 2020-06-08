/* beispiel-4.4-numerov.c   Datum:  06.05.2020 */
/* dieser Code implementiert ein vereinfachtes Numerov Problem, 
   Erweiterung notwendig bei "Korrektur notwendig fuer g!=0 !!!!" */
#pragma warning(disable : 4996)
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#ifndef M_PI
#define M_PI 3.141592653589793238
#endif

/* Routine fuer Nullstellensuche mit dem Sekantenverfahren (siehe beispiel-3.3.c)*/ 

double secant(double x1, double x2, double (*func)(double, double), int *schritt, double lambda)
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
    {         /* naechster Schatzwert x1,x2 -> xn*/
      xn=x2-func(x2,lambda)*(x2-x1)/(func(x2,lambda)-func(x1,lambda));  
      x1=x2;     /* bereite den naechsten Schritt vor:  x2 -> x1 */  
      x2=xn;     /*                                     xn -> x2 */

      (*schritt)++;      /* Schritte=Schritte+1 */
      
    }
  while(fabs(x2-x1)>tol);   /* solange Genauigkeitsziel nicht erreicht */

  return xn;   /* Gebe Nullstelle zurueck */
}


/* Code fuer die Loesung der Poisson-Gleichung mit gaussfoermiger Ladungsverteilung */

/* wegen Nullstellensuche definiere globale Felder, die in 
   Funktion zur Nullstellensuche genutzt werden koennen */

//double r0bound;                               /* definiert radius der homogenen Ladungsverteilung */
double *r_array,*g_array,*s_array,*y_loesung; /* Zeiger auf Felder fuer Numerov und Loesung */
double schrittweite_h;                        /* benutzte Schrittweite */
int num_r;                                    /* Anzahl der Stuetzstellen */ 

double gfunc(double x)
/* zu loesende Gleichung: y''(r)+g(r)*y(r)=s(r) 
   hier Definition der Funktion g(r)   */
  {
    //return 0.46332;
    //return 0.5;
    //return 0.0001713;
    //return 0.00839;   
    return 60*pow(cos(M_PI*x), 16);
  }

//Dies ist der Faktor g(x) für die Schrödinger-Gleichung mit einem zusätzlichem elektrischen Feld eps
double gfuncFeld(double x, double eps)
/* zu loesende Gleichung: y''(x)+g(x)*y(x)=0
   hier Definition der Funktion g(r)   */
{

    return 60 * pow(cos(M_PI * x), 16)+eps*x;
}

double sfunc(double x)
// zu loesende Gleichung: y''(x)+g(x)*y(x)=s(x)=0 
// hier Definition der Funktion s(r)   
{

    return 0;
}


double init_numerov(double a,double b,int n,double *r,double *g,double *s)
/* Funktionen bereitet Anwendung des Numerov Verfahrens vor 
   a,b,n sind Intervallgrenzen und Anzahl der Stuetzstellen 
   bei Verlassen der Routine beinhalten die Felder r,g und s 
   Stuetzstellen, Funktion g an den Stuetzstellen und s an den Stuetzstellen
   Rueckgabewert ist die Schrittweite h */
  {
   int i;
   double h;
   
   h=(b-a)/(n-1);                /* bestimme Intervalllaenge */
   for(i=0; i<n; i++)            /* belege die Felder mit den Werten fuer r,g und s*/
     {
       r[i]=a+i*h;
       g[i]=gfunc(r[i]);
      // printf("%lf\n", g[i]);
       s[i]=sfunc(r[i]);
     }
   
   return h;  /* Rueckgabewert ist Schrittweite h */ 
  }  

double init_numerovFeld(double a, double b, int n, double* r, double* g, double* s, double eps)
/* Funktionen bereitet Anwendung des Numerov Verfahrens vor
   a,b,n sind Intervallgrenzen und Anzahl der Stuetzstellen
   bei Verlassen der Routine beinhalten die Felder r,g und s
   Stuetzstellen, Funktion g an den Stuetzstellen und s an den Stuetzstellen
   Rueckgabewert ist die Schrittweite h, jedoch wird hier ein zusätzlicher Double Parameter eps übergeben, der die
   Stärke des elektrischen Feldes bestimmt. Dieser wird hier übergeben, nicht an gfunc sondern an gfuncFeld*/
{
    int i;
    double h;

    h = (b - a) / (n - 1);                /* bestimme Intervalllaenge */
    for (i = 0; i < n; i++)            /* belege die Felder mit den Werten fuer r,g und s*/
    {
        r[i] = a + i * h;
        g[i] = gfuncFeld(r[i], eps);
        // printf("%lf\n", g[i]);
        s[i] = sfunc(r[i]);
    }

    return h;  /* Rueckgabewert ist Schrittweite h */
}

//numerov Schritte für konstantes g(x)=lambda
void numerovdown(double *r,double *g,double *s,double h,int n,int steps,double yn,double ynm1,double *y,double lambda)
/* Funktion nutzt das Numerov Verfahren um fuer gegebene Stuetzstellen r
   und Funktionen g und s die Loesung y zu finden. 
   r, g und s sollten mit init_numerov vorbereitet werden
   h ist die Schrittweite (auch aus init_numerov)
   n ist die Anzahl der Stuetzstellen 
   steps legt fest, wieviele Numerovschritte ausgefuehrt werden
   yn und ynm1 sind die Startwerte y[n-1] und y[n-2], die Routine legt dann y[n-3] .. y[n-steps-2] fest 
   y ist ein Feld der Länge n, das die Loesungen beinhält
*/
 {
   int i;
   double fakt_u_np1,fakt_u_n,fakt_u_nm1; /* Variablen fuer Numerov Faktoren */ 
   double fakt_s;    

   y[n-1]=yn;  /* belege erste Funktionswerte mit den Startwerten  */
   y[n-2]=ynm1;
   
   for(i=n-2; i>n-steps-2; i--)   /* betrachte in den Schritten y(i-1) und y(i) um y(i+1) zu berechnen */
     {
       fakt_u_np1 = (1.0 + h * h / 12 * lambda);     /* Faktor bei y(i+1),  Korrektur notwendig fuer g!=0 !!!! */
       fakt_u_n = (1.0 - 5 * h * h / 12 * lambda);       /* Faktor bei y(i),    Korrektur notwendig fuer g!=0 !!!! */
       fakt_u_nm1 = (1.0 + h * h / 12 * lambda);     /* Faktor bei y(i-1)   Korrektur notwendig fuer g!=0 !!!! */
       fakt_s=h*h/12.0*(s[i+1]+10.0*s[i]+s[i-1]);   /* Faktor mit s */
       
       y[i-1]=(fakt_s+2.0*fakt_u_n*y[i]-fakt_u_np1*y[i+1])/fakt_u_nm1; /* Rueckwaertsiteration */   
     }

 }

//numerov Schritte für die SGL, hier wird wirklich vom Feld g gebrauch gemacht, wegen g(x)=!0
void numerovdownSGL(double* r, double* g, double* s, double h, int n, int steps, double yn, double ynm1, double* y, double lambda)
/* Funktion nutzt das Numerov Verfahren um fuer gegebene Stuetzstellen r
   und Funktionen g und s die Loesung y zu finden.
   r, g und s sollten mit init_numerov vorbereitet werden
   h ist die Schrittweite (auch aus init_numerov)
   n ist die Anzahl der Stuetzstellen
   steps legt fest, wieviele Numerovschritte ausgefuehrt werden
   yn und ynm1 sind die Startwerte y[n-1] und y[n-2], die Routine legt dann y[n-3] .. y[n-steps-2] fest
   y ist ein Feld der Länge n, das die Loesungen beinhält
*/
{
    int i;
    double fakt_u_np1, fakt_u_n, fakt_u_nm1; /* Variablen fuer Numerov Faktoren */
    double fakt_s;

    y[n - 1] = yn;  /* belege erste Funktionswerte mit den Startwerten  */
    y[n - 2] = ynm1;

    for (i = n - 2; i > n - steps - 2; i--)   /* betrachte in den Schritten y(i-1) und y(i) um y(i+1) zu berechnen */
    {
        fakt_u_np1 = (1.0 + h * h / 12 * (2*(-g[i + 1]+lambda)));     /* Faktor bei y(i+1),  Korrektur notwendig fuer g!=0 !!!! */
        fakt_u_n = (1.0 - 5 * h * h / 12 * (2*(-g[i]+lambda)));       /* Faktor bei y(i),    Korrektur notwendig fuer g!=0 !!!! */
        fakt_u_nm1 = (1.0 + h * h / 12 * (2*(-g[i -1]+lambda)));     /* Faktor bei y(i-1)   Korrektur notwendig fuer g!=0 !!!! */
        fakt_s = h * h / 12.0 * (s[i + 1] + 10.0 * s[i] + s[i - 1]);   /* Faktor mit s */

        y[i - 1] = (fakt_s + 2.0 * fakt_u_n * y[i] - fakt_u_np1 * y[i + 1]) / fakt_u_nm1; /* Rueckwaertsiteration */
    }

}
 
double f_randbed(double para,double lambda)
/* Nullstelle dieser Funktion signalisiert, dass die Loesung 
   konsistent mit der Randbedinung ist
   para ist der freie Parameter der Loesung 
   hier: para ist der Konstantewert der Loesung bei großen r 
   Die Funktion nutzt globale Variabeln, um die Felder 
   r, g, s, h und n zu erhalten   
   und feld y_loesung fuer die Loesung, zusätzlicher Double Parmater lambda wird hier benötigt*/  {      
    
    numerovdown(r_array,g_array,s_array,schrittweite_h,num_r,num_r-2,0,para,y_loesung, lambda);

    return y_loesung[0]-1;  /* zweite Randbedingung ist y[0] = 1 */
}


double f_randbedSGL(double para, double lambda)
/* Nullstelle dieser Funktion signalisiert, dass die Loesung
   konsistent mit der Randbedinung ist
   para ist der freie Parameter der Loesung
   hier: para ist der Konstantewert der Loesung bei großen r
   Die Funktion nutzt globale Variabeln, um die Felder
   r, g, s, h und n zu erhalten
   und feld y_loesung fuer die Loesung, zusätzlicher Double Parmater lambda wird hier benötigt */ {

    numerovdownSGL(r_array, g_array, s_array, schrittweite_h, num_r, num_r - 2, 0, para, y_loesung, lambda);
    return y_loesung[0];  /* zweite Randbedingung ist y[0] = 0 */
}

//Diese Funktion geht Element für Element eine Lösung von y durch und gibt den maximalen Wert zurück
double maxfinden(double *y_loesung) {
    double max;
    max = y_loesung[0];
    for (int i = 0; i < num_r; i++) {
        if (max < fabs(y_loesung[i])) {
            max = fabs(y_loesung[i]);
        }
    }
    return max;
}

//Printet eine Lösung aus
void printen(double* y) {   
    for (int i = 0; i < num_r; i++) {
        printf("%10e\n", y[i]);      
    }
    return;
}

//Rechnet den Wert Normierungsfaktor aus
double Normier(double a, double b, double *y) {
    double c = (b - a) / num_r;

    double x = 0.;
    double sum = 0.;
    for (int i = 0; i < num_r; i++) {
        sum = sum + c * (pow(y[i],2));
    }
    return 1/pow(sum,0.5);
}

int main()
  {
    double rmax,para,zero;  /* maximales r und konsistente Randbedingung und Null aus Auswertung der Randbedingungsfunktion */ 
    int i,num_schritt; /* Anzahl der Sekantenverfahren Schritte */ 
    double lambda = 0;
    int zaehl = 0;
    double* lambdas;
    int navi = 0;
    
    printf("Zu welchem Aufgabenteil wollen sie?\n1=HA 7.2-Eigenwertspektrum \n2=HA 8.1-Energieeigenwerte \n3=HA 8.2-finit size effect \n4=HA 8.3-Elektrisches Feld und Banddluecke\n");
    scanf("%d", &navi);
    num_r = 1000;
    lambdas = malloc(10 * sizeof(double));
    r_array=malloc(num_r*sizeof(double));   /* alloziere Speicher fuer r, g und s */ 
    g_array=malloc(num_r*sizeof(double));
    s_array=malloc(num_r*sizeof(double));
    y_loesung=malloc(num_r*sizeof(double));

    //schrittweite_h=init_numerov(0.0,rmax,num_r,r_array,g_array,s_array);  /* belege r,g und s mit Werten*/
    //schrittweite_h = init_numerovFeld(0.0, rmax, num_r, r_array, g_array, s_array);  /* belege r,g und s mit Werten*/


    if (navi == 1) {//-------------Aufgabe 7.2----------------
        rmax = 60;
        //FILE* data1 = fopen("xt2.txt","w");
        schrittweite_h = init_numerov(0.0, rmax, num_r, r_array, g_array, s_array);
        printf("Initialisierung beendet\n\n");
	//Schrittweise Erhöhung des Eigenwertes mit Abbruchkriterium nach dem Erreichen des 10. EW
        while ( lambda < 1 && zaehl<10) {
            //printf("Ich war hier\n");
	    //Durchführung der Shooting Methode mit Numerov 
            para = secant(-210.0, 210.0, &f_randbed, &num_schritt, lambda);
            f_randbed(para, lambda);
            //printf("%lf %lf\n", fabs(maxfinden(y_loesung) - sqrt(2)), lambda);
            //Kontrolle der maximalen Amplitude, wenn sie ->inf geht, dann ist ein EW erreicht
            if (fabs(maxfinden(y_loesung))>= 100) {
		    //Kontrolle ob es zwei verschiedene EW sind
                if(lambdas[zaehl-1]+0.00033<lambda ||zaehl==0){
		    //EW ausgeben	
                    printf("Hab einen\n");

                    printf("%10e\n", lambda);

                    lambdas[zaehl] = lambda;

                    zaehl++;
                }
                
                
            }
            lambda = lambda + 0.00001;
        }
        //fclose(data1);
	//Finale Aufzählung der EW
        for (i = 0; i < 10; i++) {
            printf("Eigenwert %d = %lf\n", i + 1, lambdas[i]);
        }
    }
    else if (navi == 2) {//-------------Aufgabe 8.1----------------
        rmax = 8;
	//Initialisierung
        schrittweite_h = init_numerov(0.0, rmax, num_r, r_array, g_array, s_array);
        printf("Initialisierung beendet\n\n");
        //FILE* data = fopen("xt.txt", "w");
	//Zählvariable g für die NR des EW 
        int g = 1;
        for (double i = 0; i < 60; i = i + 0.01) {
            //para = secant(0.1, 1, &f_randbedSGL, &num_schritt, i);
	    //Ausführung des Numerov Verfahrens mit Wert i=lambda=EW	
            f_randbedSGL(0.000001, i);
	    //Ist die Randbedingung erfüllt?
            if (fabs(y_loesung[0]) <= 10e-7) {
		//Ausgabe eines gefundenen EWs und der dazugehörigen Normierungskonstante
                printf("Energiewert %d: %lf\n", g, i);
                printf("Normierungskonstante N=%10e\n", Normier(0, 8, y_loesung));
                //fprintf(data, "%d %lf\n", g, i);
                g++;
            }

        }
        //f_randbedSGL(0.000001, 19.660000);
        //printen(y_loesung);
        //fclose(data);
        
    }

    else if (navi == 3) {//-------------Aufgabe 8.2----------------
	//Annahme des Parameters L um den f-s-e zu untersuchen
        printf("Wir brauchen hier eine neue Länge:");
        scanf("%lf", &rmax);
	//Numerov Initialisierung
        schrittweite_h = init_numerov(0.0, rmax, num_r, r_array, g_array, s_array);

        printf("Initialisierung beendet\n\n");
	//Zählvariable g für die NR des EW 
        int g = 1;
        for (double i = 0; i < 30 && g<=25; i = i + 0.01) {
            //para = secant(0.1, 1, &f_randbedSGL, &num_schritt, i);
	    //Ausführung des Numerov Verfahrens mit Wert i=lambda=EW
            f_randbedSGL(0.000001, i);
            //Ist die Randbedingung erfüllt?		
            if(fabs(y_loesung[0]) <= 10e-7) {
		//Ausgabe eines gefundenen EWs und der dazugehörigen Normierungskonstante
                printf("Energiewert %d: %lf\n", g, i);
                g++;
            }

        }



    }
    else if (navi == 4) {//-------------Aufgabe 8.3----------------
        rmax = 8;

        //FILE* data = fopen("xt2.txt", "w");
        int g = 1;//EW Zähler
        int abb = 0;//Bandlücken Zähler
        double zwi = 0;//Vorheriger EW wird hier abgespeichert
	//Hier gilt es das elektrische Feld zu erhöhen
        for (double eps = 0; eps < 0.3; eps = eps + 0.01) {
		
            abb = 0;//Bandlücken Zähler auf 0 setzen
            //Neu Initialisierung von Numerov für ein elektrisches Feld
            schrittweite_h = init_numerovFeld(0.0, rmax, num_r, r_array, g_array, s_array, eps);
	    
            printf("\n\nInitialisierung fuer E-Feld=%lf beendet\n",eps);
            //fprintf(data, "%lf", eps);
            zwi = 0;
	    //Durchführung des vorherigen Algorithmuse, also Schritt für Schritt den EW erhöhen	
            for (double i = 0; i < 60 && abb < 2; i = i + 0.001) {
                f_randbedSGL(0.000001, i);
                if (fabs(y_loesung[0]) <= 10e-6) {
                    //printf("%lf\n", i);
		    //Hier wird jetzt auf eine Bandlücke überprüft indem der Wert des vorherigen und des neuen EW verglichen wird 		
                    if (i - zwi >  3 && abb < 2 && zwi != 0) {
			//Werte der ersten zwei Bandlücke werden notiert    
                        printf("Energiewerte %lf: %lf, Diff %lf\n", zwi, i, fabs(i - zwi));
                        //fprintf(data, " %lf", fabs(i - zwi));
			//Bandlücken-Zähler erhöhen
                        abb = abb + 1;
                    }
                    zwi = i;
                    g++;
                }

            }
            g = 1;
            //fprintf(data, "\n", eps);
        }
        //fclose(data);
    }

    free(lambdas);
    free(r_array);
    free(g_array);
    free(s_array);
    free(y_loesung);

    return 0;
  }
