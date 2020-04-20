#include<stdio.h>
#include<stdlib.h>
void printarray(int *array,int length){
  int i;
  for (i=0;i<length;++i)
    printf("%d; ",array[i]);
  printf("\n");
}
void initializieren(int *array,int length,int wert){
  int i;
  for (i=0;i<length;++i)
    array[i]=wert;
}
void rotieren1(int *array,int length){
  int i, temp;
  temp = array[0];
  for (i = 0; i < length-1; i++)
     array[i] = array[i+1];
  array[i] = temp;
}
void rotierenk(int *array,int length,int k){
  int i;
  for (i=0;i<k;++i)
    rotieren1(array,length);
}
void umdrehen(int *array, int length)
{
  int temp,start,end;
  for (start=0,end=length-1;start<end;start++,end--){
    temp = array[start];   
    array[start] = array[end];
    array[end] = temp;
  }
}
int suchen(int *a,int length,int *b,int lengthb){
  int ret=0;
  int i,j;
  //Für jeden a[i] versuchen wir a[i]=b[0] überprüfen
  for (i=0;i<=length-lengthb;++i){
    for (j=0;j<lengthb;++j)
    //Wenn es nicht klappt machen wir einen Schritt weiter
      if (a[i+j]!=b[j])
       break;
    if (j==lengthb)
    //Wenn es klappt geben wir die Pozition zurück
      return i;
  }
  return -1; 
} 
int main(int argc,char *argv[]){
  int array[]={1,2,3,4};
  int A[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  int B[3]= {4,5,6};
  int C[2]={5,7};
  int D[2]={9,10};
  int ergebnis;
  printarray(array,4);
  rotierenk(array,4,2);
  printarray(array,4);
  umdrehen(array,4);
  printarray(array,4);
  ergebnis=suchen(A,10,B,3);
  printf("A in B ist %d\n", ergebnis);
  ergebnis=suchen(A,10,C,2);
  printf("A in C ist %d\n", ergebnis);
  ergebnis=suchen(A,10,D,2);
  printf("A in D ist %d\n", ergebnis);
}
