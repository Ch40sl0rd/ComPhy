#include "AblT.h"
double kinEnerg( double *x, double m, unsigned int t, unsigned int tmax ){
  double v = AblT(x, t, tmax);
  return 0.5*m*v*v;
}
