#include<stdio.h>
#include<math.h>
int sgn_own(double x){
  if (x==0)
   return 0;
  if (x>0)
   return +1;
  else
   return -1;
}
double betrag_own( double x ){
  if (x>=0)
    return x;
  else
    return -x;
}
double cos_own( double x ){
  double ret=0;
  int i;
  double ai=1;
  double aip1;
  //cosx=\sum_i a_i
  //a_0=1.
  //a_i+1=(-1)*x**2/((2*i+1)*(2*i+2))
  for (i=0;i<100;++i){
    ret+=ai;
    aip1=ai*(-1.)*x*x/((2.*i+1.)*(2.*i+2.));
    ai=aip1;
  }
  return ret;
}
double wurzel_own(double c){
  double x;
  double y;
  x=5.;//Anfangs wert
  y=0.5*(x+c/x);//Rekurzion Anfangen
  do{
   x=0.5*(y+c/y);
   if (betrag_own(x-y)<0.000005)//Wenn wir den Fixpunkt erreichen, sind wir Fertig
    break;
   //Sonst wiederholen wir
   y=x;
  }while(1);
  return x;
}
int main(int argc, char *argv[]){
  double c=wurzel_own(0.5);
  printf("%le\n",c);
}
