#include<stdio.h>
#include<stdlib.h>
double fown(double x){
  return x;
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
int main(){
  funktiondodatei(0.01,1.,2.,&fown,"temporary");
}

