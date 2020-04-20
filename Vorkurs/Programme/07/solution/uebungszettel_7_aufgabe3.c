#include<stdio.h>
#include<stdlib.h>
#include<math.h>
double fown(double x){
  return x*x;
}
double funktiondodatei(double s,double x1, double x2,double (*f)(),char *filename){
  FILE *out=fopen(filename,"w");
  if (out== NULL){
    printf("%s file kann nicht eroffnet werden\n",filename);
    exit(1);
  }
  int i;
  int nbins;
  nbins=(int)(x2-x1)/s;
  for (i=0;i<nbins;++i){
    double temp=f(x1+i*s);
    fprintf(out,"%e\t%e\n", x1+i*s,f(x1+i*s));
  }
  fclose(out);
}
double funktionintegratea(double a,double b,double (*f)(double),unsigned int n){

  double summe=0;
  double deltax=(b-a)/(double)n;
  int i;
  for (i=0;i<n;++i){
    summe+=(f(i*deltax)+f((i+1.)*deltax))/2.;
  }
  summe*=deltax;
}

double funktionintegrateb(double a,double b,double (*f)(double),double e){

  double summeold=0;
  double summenew=1;
  int n;
  for (n=10;fabs(summenew-summeold)/fabs(summenew)>e;n*=2){
    summeold=summenew;
    summenew=0.0;
    double deltax=(b-a)/(double)n;
    int i;
    for (i=0;i<n;++i){
      summenew+=(f(i*deltax)+f((i+1.)*deltax));
    }
    summenew*=deltax;
    summenew/=2.;
  }
  return summenew;
}

int main(){
  double sum=funktionintegrateb(0,1,&fown,0.0001);
  printf("Integral of f ist %e\n",sum);
//  funktiondodatei(0.01,1.,2.,&fown,"temporary");
}

