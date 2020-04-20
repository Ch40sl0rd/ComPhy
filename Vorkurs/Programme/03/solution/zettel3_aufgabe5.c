#include<stdio.h>
#include<stdlib.h>
#include<math.h>
double riemannzeta_own( double x ){
  double ret=0;
  int i;
  double ai=1;
  double aip1;
  //xis=\sum_i a_i
  //a_1=1.
  //a_i+1=1/(1+1/i)**x*ai
  for (i=1;i<1000;++i){
    ret+=ai;
    aip1=ai*pow(1./(1.+1./i),x);
    ai=aip1;
  }
  return ret;
}
int main(int argc,char *argv[]){
  double x=riemannzeta_own( 2. );
  printf("Riemannzeta %e\n", x);
}
