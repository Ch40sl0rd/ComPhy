/***********************************************************************
 *
 * compile with
 *
 *   gcc -O3 -Wall -pedantic A2_trapez.c -o A2_trapez -lm
 *
 * run with
 *   ./A2_trapez
 *
 ***********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>  /* needed for memcpy */
#include <math.h>

/* useful, stringification of macro */
#define str(_a) #_a
#define xstr(_a) str(_a)

/* definiere per default to 8 , also double precision*/
#define PRECISION_SIZE 8
#define PRECISION_TYPE double

/* function pointer type definition;
 * f vom Typ func_type mit Argument vom Typ PRECISION_TYPE, Parameterliste p 
 * Return-Typ PRECISION_TYPE 
 *
 * Dieser Funktionstyp wird nachher als Integrand verwendet
 */
typedef PRECISION_TYPE (*func_type ) ( PRECISION_TYPE const x, void * p );

/***********************************************************************
 * elementarer Trapze-Schritt
 *
 * Brauchen wir nachher nicht explizit
 ***********************************************************************/
PRECISION_TYPE trapez ( func_type f, PRECISION_TYPE const xa, PRECISION_TYPE const xe, void * p ) {

  return( ( f(xe, p) + f(xa, p) ) * 0.5 * ( xe - xa ) );

} /* end of trapez */

/***********************************************************************
 * Trapez-Integration
 *
 * nur Beitraege der Zwischenwerte xa, xa+h, ..., xe
 * Parameterliste 
 ***********************************************************************/
PRECISION_TYPE trapez_integration_restarted ( func_type f, PRECISION_TYPE const xa, PRECISION_TYPE const xe, PRECISION_TYPE const h, void * p ) {

  /* Anzahl der Zwischenschritte */
  int const nstep = ( xe - xa ) / h;

  /* fprintf( stdout, "# [trapez_integration_restarted] nstep = %4d\n", nstep ); */

  /* Initialisierung */
  PRECISION_TYPE tsum = 0.;
  PRECISION_TYPE x    = xa;

  /* Beitraege der Zwischenwerte */
  for ( int i = 0; i <= nstep ; i++ ) {
    /* fprintf ( stdout, "# [trapez_integratio_restarted] h = %e x = %e f = %e\n", h, x, f(x,p) ); */
    tsum += f( x, p );
    x += h;
  }
  /* mal Schrittweite */
  tsum *= h;

  return( tsum );

} /* end of trapez_integration */

/***********************************************************************
 * Beispielfunktionen zum Integrieren
 ***********************************************************************/

/***********************************************************************
 * Bsp. Cauchy-Verteilungsfunktion
 *
 * s / ( s^2 + / ( x - t )^2 ) / pi
 ***********************************************************************/
PRECISION_TYPE cauchy ( PRECISION_TYPE x , void * p ) {
  PRECISION_TYPE const s = ((PRECISION_TYPE*)p)[0];
  PRECISION_TYPE const t = ((PRECISION_TYPE*)p)[1];
  return ( s / ( s * s + ( x - t ) * ( x - t ) ) / M_PI  );
}  /* end of cauchy */

/***********************************************************************
 * Bsp. Parabel-Funktion
 *
 * a x^2 + b x + c
 ***********************************************************************/
PRECISION_TYPE parabola ( PRECISION_TYPE x , void * p ) {
  PRECISION_TYPE const a = ((PRECISION_TYPE*)p)[0];
  PRECISION_TYPE const b = ((PRECISION_TYPE*)p)[1];
  PRECISION_TYPE const c = ((PRECISION_TYPE*)p)[2];
  return ( ( a * x + b ) * x +c );
}  /* end of parabola */

/***********************************************************************
 * Bsp. Integrand fuer Gamma-Funktion
 *
 * exp ( -t ) t^( z - 1 ), z reell
 ***********************************************************************/
PRECISION_TYPE Gamma_integrand ( PRECISION_TYPE x, void * p  ) {
  PRECISION_TYPE const z = ((PRECISION_TYPE*)p)[0];
  return( 
      /* Funktionsvariante in Abh. vom Genauigkeitstyp */
#if   PRECISION_SIZE == 4
      /* single precison / float  */
      expf ( -x ) * powf ( x, z - 1. )
#elif PRECISION_SIZE == 8
      /* doubpe precision / double */
      exp  ( -x ) * pow  ( x, z - 1. )
#endif
  );
}  /* end of Gamma_integrand */

/***********************************************************************
 * MAIN PROGRAM
 ***********************************************************************/
int main(int argc, char **argv) {

  /* Integrand und Parameter fuer Integration von Gamma_integrand
   *
   * brauchen nur fuer Parameter 1 <= p <= 2 integrieren, sonst
   * p Gamma ( p ) = Gamma ( p + 1 ) ausnutzen
   *
   */ 
  func_type f = Gamma_integrand;
  PRECISION_TYPE p[1] = { 3.5 };       /* Argument z fuer Gamma ( z ) */
  PRECISION_TYPE int_xa = 0.;          /* untere Intervallgrenze */
  PRECISION_TYPE int_xe = 20.;         /* obere Intervallgrenze  */

  PRECISION_TYPE int_epsabs = 1.e-08;  /* geforderete max. absolute Aenderung */
  PRECISION_TYPE int_epsrel = 1.e-06;  /* geforderete max. relative Aenderung */

  unsigned int int_nstep_min = 10;     /* minimale Anzahl Schritte, werden auf jeden Fall durchgefuehrt */
  unsigned int int_nstep_max = 40;     /* maximal Anzahl Schritte, abh. von der Konvergenz */

  /* Startwerte */
  PRECISION_TYPE diffrel = 2 * int_epsrel;
  PRECISION_TYPE diffabs = 2 * int_epsabs;

  PRECISION_TYPE h = int_xe - int_xa;

  unsigned int nstep = 1;
  
  /* Start der Integrations-Iteration mit elementarem Trapez-Schritt */
  PRECISION_TYPE int_val = trapez ( f, int_xa, int_xe, p );  /* nur Beitraege xa, xe */

  /* fprintf ( stdout, "# [main] nstep = %3u h = %25.16e   int_val = %25.16e\n", nstep, int_xe - int_xa, int_val ); */
  nstep++;

    /* Iteration bis hinauf zu nstep_min */
  for ( ; nstep <= int_nstep_min; ) {
    PRECISION_TYPE const xa = int_xa + h / 2.;
    PRECISION_TYPE const xe = int_xe - h / 2.;
    int_val = 0.5 * ( int_val + trapez_integration_restarted ( f, xa, xe, h, p ) );
    h *= 0.5;
    /* fprintf ( stdout, "# [main] nstep = %3u   xa = %e xe = %e   h = %25.16e   int_val = %25.16e\n", nstep, h, int_val ); */
    nstep++;
  }
  
  /*
   * Iteration bis zu max. nstep_max
   *
   * Stop-Kriterum: nstep_max, epsrel, epsabs
   *
   */
  while ( nstep <= int_nstep_max && ( diffrel > int_epsrel || diffabs > int_epsabs ) ) {

    PRECISION_TYPE const xa = int_xa + h / 2.;
    PRECISION_TYPE const xe = int_xe - h / 2.; 
    PRECISION_TYPE val_new = 0.5 * ( int_val + trapez_integration_restarted ( f, xa, xe, h, p ) );
    h *= 0.5;

    diffabs = fabs( int_val - val_new );
    diffrel = diffabs / fabs( int_val + val_new ) * 2.;
    
    /* aktualisiere int_val */
    int_val = val_new;

    /* fprintf ( stdout, "# [main] nstep = %3u h = %25.16e   int_val = %25.16e   eps %e %e\n", nstep, h, int_val, diffabs, diffrel ); */
    
    nstep++;

  }  /* end of nstep, eps iteration */
    
  nstep--;
  fprintf ( stdout, "# [main] nstep = %3u h = %25.16e   int_val = %25.16e   eps %e %e\n", nstep, h, int_val, diffabs, diffrel );

  return ( 0 );
}
