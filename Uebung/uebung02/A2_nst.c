/***********************************************************************
 *
 * compile with
 *
 *   gcc -Wall -pedantic A2_nst.c -o A2_nst -lm
 *
 * run with
 *   ./A2_nst
 *
 ***********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* useful, stringification */
#define str(_a) #_a
#define xstr(_a) str(_a)

/* if PRECISION_SIZE is not defined during integration, define per default to 8 */
#define PRECISION_SIZE 8
#define PRECISION_TYPE double

/* function point type definition;
 * f or type func_type will take an argument x and return a value, both of type PRECISION_TYPE */
typedef PRECISION_TYPE (*func_type ) ( PRECISION_TYPE const x);

/***********************************************************************
 * first, symmetric derivative
 ***********************************************************************/
PRECISION_TYPE deriv_sym_1 ( func_type f, PRECISION_TYPE const x, PRECISION_TYPE const h ) {

  return( ( f(x+h) - f(x-h) ) / (2.*h) );

} /* end of deriv_sym_1 */

/***********************************************************************
 * second, symmetric derivative
 ***********************************************************************/
PRECISION_TYPE deriv_sym_2 ( func_type f, PRECISION_TYPE const x, PRECISION_TYPE const h ) {

  return( ( f(x+h) + f(x-h) - 2 * f(x) ) / (h*h) );

} /* end of deriv_sym_2 */

PRECISION_TYPE myfunc ( PRECISION_TYPE x ) {
  return ( pow( log( x ), 1./x ) );
}

PRECISION_TYPE secant_deriv (func_type f, PRECISION_TYPE const x1, PRECISION_TYPE const x2, PRECISION_TYPE const h ) {

  PRECISION_TYPE const df1 = deriv_sym_1 ( f, x1, h );
  PRECISION_TYPE const df2 = deriv_sym_1 ( f, x2, h );

  return( x1  - df1 * ( x1 - x2 ) / ( df1 - df2 ) );
}


/***********************************************************************
 * MAIN PROGRAM
 ***********************************************************************/

int main(int argc, char **argv) {

  /* */ 
  func_type f = myfunc;


  /* parameters for secant method */
  PRECISION_TYPE nst_x1  =  4.;
  PRECISION_TYPE nst_x2  =  8.;
  PRECISION_TYPE nst_epsabs = 1.e-08;
  PRECISION_TYPE nst_epsrel = 1.e-06;
  unsigned int nst_maxiter = 10000;

  PRECISION_TYPE deriv_h0 = 1.0;
  unsigned int deriv_nstep = 32;
  PRECISION_TYPE deriv_dstep = 0.5;

  /* output file for iteration on stepsize h */
  char filename[400];
  sprintf(filename, "nst.%s", xstr(PRECISION_TYPE) ) ;
  FILE * ofs = fopen( filename, "w" );

  /* interation on step sizes for function evaluation around x0 */
  PRECISION_TYPE h = deriv_h0;

  for ( unsigned int i = 0; i <= deriv_nstep; i++ ) {

    PRECISION_TYPE x1 = nst_x1;
    PRECISION_TYPE x2 = nst_x2;

    /* variables for absolute and relative differenc */
    PRECISION_TYPE dabs = fabs( x1 - x2 );
    PRECISION_TYPE drel = fabs( x1 - x2 ) / fabs( x1 + x2 ) * 2.;
    unsigned int count = 0;

    /* loop on condition on rel./abs. change in iteration and iteration count below max. allowed number */
    while ( ( dabs > nst_epsabs || drel > nst_epsrel ) && count < nst_maxiter ) {

      /*update x1, x2 ; via tmp */
      PRECISION_TYPE tmp  = secant_deriv ( f, x1, x2 , h );
      x2 = x1;
      x1 = tmp;

      /* new absolute and relative difference */
      dabs = fabs( x1 - x2 );
      drel = fabs( x1 - x2 ) / fabs( x1 + x2 ) * 2.;

      /* increase the iteration counter */
      count++;

      /* show result during secant method iteration */
      fprintf ( stdout, "%25.16e   %25.16e %25.16e    %25.16e   %25.16e   %8u\n", h, x1, x2, f ( x1 ), deriv_sym_1 ( f, x1, h ), count );
    }

    /* save results to file */
    fprintf ( ofs, "%25.16e   %25.16e %25.16e    %25.16e   %25.16e   %8u\n", h, x1, x2, f( x1 ), deriv_sym_1 ( f, x1, h ),  count );
    fprintf ( ofs, "# =================================================\n" );

    /* decrease h by factor deriv_dstep */
    h *= deriv_dstep;
  }

  /* close file */
  fclose ( ofs );

  return ( 0 );
}
