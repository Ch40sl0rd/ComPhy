#include<stdio.h>
#include<stdlib.h>
#include<math.h>
double square_to( double *a){
  double temp=(*a)*(*a);
  (*a)=temp;
  return temp;
}
double root_to( double *a){
  double temp=sqrt((*a));
  (*a)=temp;
  return temp;
}

int main(int argc,char *argv[]){
  double m=10;
  double x=square_to(&m);
  printf("m=%le n=%le\n", x,m);
  x=root_to(&m);
  printf("m=%le n=%le\n", x,m);
}        
