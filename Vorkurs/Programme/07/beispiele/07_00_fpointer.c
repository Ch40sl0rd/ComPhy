#include <math.h>
#include <stdio.h>

int main(void){

  double (*fp)(double);

  fp = cos;

  printf("%f\n", (*fp)(0.3) );

  fp = sin;

  printf("%f\n", (*fp)(0.3) );

  return 0;
}
