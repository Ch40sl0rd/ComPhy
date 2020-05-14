void RuKu_3 ( int nDifEqu, /* # der Differentialgleichungen */
              double h, /* Schrittweite */
              double t, /* Kurvenparameter */
              double y[], /* Bahnkurve [nDifEqu] */  
              double yh[], double k1[], double k2[], double k3[], /* Hilfsfelder [nDifEqu] */
              void (*derhs) ( int, double, double[], double[] ) /* Funktion zur Berechnung der rechten Seite */
    ){

/* Berechnet 1 Runge-Kutta Schritt (3. Ordnung) zur Loesung des DG-Systems:
   y'(t) = f(y(t),t); y \in R^n */

/* Variablen */
  double h2;
  int i;

  h2 = 0.5 * h;
  (*derhs)( nDifEqu, t   , y , k1 );
  for (i=0; i<nDifEqu; i++) { yh[i] = y[i] + h2 * k1[i]; }
  (*derhs)( nDifEqu, t+h2, yh, k2 );
  for (i=0; i<nDifEqu; i++) { yh[i] = y[i] - h  * k1[i] + 2 * h * k2[i]; }
  (*derhs)( nDifEqu, t+h2, yh, k3 );
  for (i=0; i<nDifEqu; i++) { y[i] += 
          ( h2 * (k1[i]+k3[i]) + 2 * h * k2[i] ) / 3; }

}



