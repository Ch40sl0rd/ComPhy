#include <stdlib.h>
#include <stdio.h>
double AblT( double* x, unsigned int t, unsigned int tmax ){
  if( t+1 >= tmax ){
    printf("x hat nur %d Elemente!\n", tmax);
    exit(55);
  }
  double dt = 0.01;
  return (x[t+1]-x[t])/dt;
}
