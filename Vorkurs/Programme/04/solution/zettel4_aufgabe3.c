#include<stdio.h>
#include<stdlib.h>
#include<math.h>
void solve(double a,double b,double c,double *x1,double *x2){
  if ((b*b-4+a*c)<0){
    printf("Es gibt keine Lösungen\n");
  }
  if ((b*b-4*a*c)==0.){
    *x1=-b/(2.*a);
    *x2=-b/(2.*a);
  }
  else{
    *x1=(-b+sqrt(b*b-4*a*c))/(2.*a);
    *x2=(-b+sqrt(b*b-4*a*c))/(2.*a);
  }
}

int main(int argc,char *argv[]){
  double a=1.,b=4.,c=4.;
  double x1,x2;
  solve(a,b,c,&x1,&x2);
  printf("x1=%le x2=%le\n", x1,x2);
}        
