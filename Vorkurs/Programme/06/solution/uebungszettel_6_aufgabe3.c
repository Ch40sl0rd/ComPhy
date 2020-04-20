#include <stdio.h>
#include <stdlib.h>
int *merge(int *a, int n, int *b, int m){
   int *newlist=(int*)malloc(sizeof(int)*(n+m));
   if (newlist ==NULL){
     printf("Wir haben leider nicht genug Speicher\n");
     exit(1);
   }
   int index1,index2;
   int i;
   for (i=0,index1=0,index2=0;(index1<n) && (index2<m);++i){
     if (a[index1] < b[index2]){
      newlist[i]=a[index1];
      index1++;
     }
     else{
      newlist[i]=b[index2];
      index2++;
     }
   }
   if ((index1 == n) && (index2 <m)){
     for(;index2<m;++index2,++i)
       newlist[i]=b[index2];
   }
   if ((index2 == m) && (index1 <n)){
     for(;index1<n;++index1,++i)
       newlist[i]=a[index1];
   }
   return newlist;
}
void mergesort(int *a, int n){

   int *newlist;
   int i=0;
   if (n==0 || n==1)
     return;
   //in array in zwei Teilarrays zerlegen
   mergesort(a,n/2);
   mergesort(&a[n/2],n-n/2);
   //Beide Teilarrays sind schon sortiert
   //Wir mussen sie nur kombinieren
   newlist=merge(a,n/2,&a[n/2],n-n/2);
   for (i=0;i<n;++i)
     a[i]=newlist[i];
   free(newlist);
}
int main() { 
   int i;
   int *a=(int*)malloc(sizeof(int)*10);
   if (a==NULL){
     printf("Wir haben leider nicht genug Speicher\n");
     exit(1);
   }
   a[0]=10;
   a[1]=14;
   a[2]=19;
   a[3]=26;
   a[4]=27,
   a[5]=31;
   a[6]=33;
   a[7]=42;
   a[8]=44;
   a[9]=0;
   printf("List before sorting\n");
   for(i = 0; i <10; i++)
      printf("%d ", a[i]);
    fflush(stdout);
   mergesort(a, 10);
   printf("\nList after sorting\n");
   for(i = 0; i <10; i++)
      printf("%d ", a[i]);
   printf("\n");
   free(a);
}
