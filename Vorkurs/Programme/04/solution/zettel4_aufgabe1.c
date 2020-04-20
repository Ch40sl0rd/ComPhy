#include<stdio.h>
#include<stdlib.h>
void tauschen( int *a, int *b){
  int temp=(*a);
  *a=*b;
  *b=temp;
}
int main(int argc,char *argv[]){
  int m=10;
  int n=20;
  printf("m=%d n=%d\n", m,n);
  tauschen(&m,&n);
  printf("m=%d n=%d\n", m,n);
}
