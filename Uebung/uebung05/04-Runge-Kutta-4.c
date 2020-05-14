void RuKu_4 ( int nDifEqu, /* # der Differentialgleichungen */
              double h,    /* Schrittweite */
              double t,    /* Kurvenparameter */
              double y[],  /* Bahnkurve [nDifEqu] */ 
                           /* Eingabe: y(t) */
                           /* Rueckgabe; y(t+h) */ 
              double yh[], /* Hilfsfeld [nDifEqu] */
              /* Hilfsfelder [nDifEqu]: */
              double k1[], double k2[], double k3[], double k4[],
              /* Funktion zur Berechnung der rechten Seite: */ 
              void (*derhs) ( int, double, double[], double[] ) 
    ){

/* Berechnet 1 Runge-Kutta Schritt (4. Ordnung) zur Loesung des DG-Systems:
   y'(t) = f(y(t),t); y \in R^n */

/* Variablen */
  double h2;
  int i;

  h2 = 0.5 * h;
  (*derhs)( nDifEqu, t   , y , k1 );
  for (i=0; i<nDifEqu; i++) { 
      yh[i] = y[i] + h2 * k1[i]; 
  }
  (*derhs)( nDifEqu, t+h2, yh, k2 );
  for (i=0; i<nDifEqu; i++) { 
      yh[i] = y[i] + h2 * k2[i]; 
  }
  (*derhs)( nDifEqu, t+h2, yh, k3 );
  for (i=0; i<nDifEqu; i++) { 
      yh[i] = y[i] + h  * k3[i]; 
  }
  (*derhs)( nDifEqu, t+h , yh, k4 );
  for (i=0; i<nDifEqu; i++) { 
      y[i] += ( h2 * (k1[i]+k4[i]) + h * (k2[i]+k3[i]) ) / 3; 
  }

}
