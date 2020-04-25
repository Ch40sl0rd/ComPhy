/***********************************************************************
 *
 * compile with
 *
 *   gcc -DPRECISION_SIZE=<4/8>  -Wall -pedantic A2_abl.c -o A2_abl -lm
 *     for single/double precision
 *
 *   gcc -DPRECISION_SIZE=16  -Wall -pedantic A2_abl.c -o A2_abl -lquadmath -lm
 *     for quadruple precision ( using libquadmath )
 *
 * run with
 *   ./A2_abl
 *
 ***********************************************************************/
#if PRECISION_SIZE == 16
#include <quadmath.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* useful, stringification */
#define str(_a) #_a
#define xstr(_a) str(_a)

/* if PRECISION_SIZE is not defined during integration, define per default to 8 */
#ifndef PRECISION_SIZE
#define PRECISION_SIZE 8
#endif

/* depending on PRECISION_SIZE, define the data type and format for formatted output */
#if   PRECISION_SIZE ==  4
#  define PRECISION_TYPE float
#  define FORMAT_STR "%25.16e"
#elif PRECISION_SIZE ==  8
#  define PRECISION_TYPE double
#  define FORMAT_STR "%25.16e"
#elif PRECISION_SIZE == 16
#  define PRECISION_TYPE __float128
#  define FORMAT_STR "%25.16Qe"
#else
#  error "unrecognized precision byte size"
#endif

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

/***********************************************************************
 * third, symmetric derivative
 ***********************************************************************/
PRECISION_TYPE deriv_sym_3 ( func_type f, PRECISION_TYPE const x, PRECISION_TYPE const h ) {

  return( ( ( f(x+2*h) - f(x-2*h) ) - 2 * ( f(x+h) - f(x-h) ) )  / (2*h*h*h) );

} /* end of deriv_sym_3 */


/***********************************************************************
 * find eps, such that 1 + eps == 1 to machine precision
 * with n bit mantissa should be 1/2^(n+1)
 ***********************************************************************/
PRECISION_TYPE find_eps ( void ) {

  PRECISION_TYPE const one  = 1.;    /* const value 1 in data type */
  PRECISION_TYPE const decr = 0.5;   /* decrement factor 1/2       */
  PRECISION_TYPE eps = 1.;           /* start value for interation */

  while ( (one + eps) - one > 0 ) {
    eps *= decr;
  }
  /* format string for eps according to data type */
  fprintf ( stdout, "# [find_eps] %s eps = "FORMAT_STR"\n", xstr(PRECISION_TYPE), eps );

  return( eps );
}


/***********************************************************************
 * MAIN PROGRAM
 ***********************************************************************/

int main(int argc, char **argv) {

/* select logarithm function from math ( / quadmath ) library with appropriate argument and return type */
#if   PRECISION_SIZE ==  4
  func_type f = logf;
#elif PRECISION_SIZE ==  8
  func_type f = log;
#elif PRECISION_SIZE == 16
  func_type f = logq;
#endif

  /* call to find_eps */
  find_eps ();

  /* 1st/2nd/3rd derivative of log function at x0 = 3 */
  PRECISION_TYPE x0 = 3.;
  PRECISION_TYPE h0 = 1.;
  unsigned int n_step = 1000;
  PRECISION_TYPE dec_step = 0.95;

  /* output file for iteration on stepsize h */
  char filename[400];
  sprintf(filename, "deriv_sym.%s", xstr(PRECISION_TYPE) ) ;
  FILE * ofs = fopen( filename, "w" );

  /* interation on step sizes for function evaluation around x0 */
  PRECISION_TYPE h = h0;
  for ( unsigned int i = 0; i <= n_step; i++ ) {

    /* evaluate function f ( = log ) and 1st,2nd,3rd derivative */
    PRECISION_TYPE fderiv[4] = {
      f(x0),
      deriv_sym_1 ( f, x0, h ),
      deriv_sym_2 ( f, x0, h ),
      deriv_sym_3 ( f, x0, h ) 
    };

    /* save results to file */
    fprintf ( ofs, FORMAT_STR"    "FORMAT_STR"    "FORMAT_STR" "FORMAT_STR" "FORMAT_STR" "FORMAT_STR"\n", x0, h, 
       fderiv[0], fderiv[1], fderiv[2], fderiv[3] );

    /* decrease h by factor dec_step */
    h *= dec_step;
  }

  /* close file */
  fclose ( ofs );

  return ( 0 );
}
