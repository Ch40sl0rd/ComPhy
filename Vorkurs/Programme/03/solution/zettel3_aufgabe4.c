#include<stdio.h>
#include<stdlib.h>
#include<math.h>
double exp_own( double x ){
  double ret=0;
  int i;
  double ai=1;
  double aip1;
  //exp x=\sum_i a_i
  //a_0=1.
  //a_i+1=(+1)*x/((i+1))
  for (i=0;ai>0.00000001;++i){
    ret+=ai;
    aip1=ai*(+1.)*x/((i+1.));
    ai=aip1;
  }
  return ret;
}
double ln_own( double x ){
  double ret=0;
  int i;
  double ai=2.*(x-1.)/(x+1.);
  double aip1;
  //ln x=2*\sum_i a_i
  //a_0=(x-1)/(x+1)
  //a_i+1=((x-1)/(x+1))**2*(2*i+1)/(2*i+3):
  for (i=0;i<1000;++i){
    ret+=ai;
    aip1=ai*pow((x-1.)/(x+1.),2.*(2*i+1))/(2*i+3);
    ai=aip1;
  }
  return ret;
}
double power(double x,double y){
  double z=ln_own(x);
  return exp_own(y*z);
}
int main(int argc,char *argv[]){
  double x,y;
  x=exp_own(1.);
  y=ln_own(x);
  printf("y=%e\n",ln_own(2.714));
  y=power(3,4.);
  printf("y=%e\n",y);
}
