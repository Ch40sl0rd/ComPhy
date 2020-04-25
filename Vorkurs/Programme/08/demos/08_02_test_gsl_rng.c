#include <stdio.h>
#include <gsl/gsl_rng.h>
int main (void) {
  gsl_rng_type const * generator_type;
  gsl_rng * r;
  gsl_rng_env_setup();
  generator_type = gsl_rng_ranlxd2;
  r = gsl_rng_alloc(generator_type);
  gsl_rng_set(r, 12345);
  double u = 0.0;
  for(int i = 0; i < 50; i++) {
      u = gsl_rng_uniform (r);
      printf ("%.5f\n", u);
  }
  gsl_rng_free(r);
  return 0;
}
