#include<stdio.h>
#include"arrayhelpers.h"

int main(){
  double *array;
  int n,i;

  n=42;
  /*Hole Speicher fuer n Eintrage*/
  array=(double*)malloc(sizeof(double)*n); //Fehler1: malloc(n) alokkiert nur "n" bytes,
                                           //Wir breauchen n*8 bytes
  if (array==NULL){
    printf("Fehler, nicht genug Speicher\n");
    exit(1); //Fehler2: Wenn nicht genug Spiecher gibt, wir müssen sofort abbrechen
  }
  for (i=0;i<n;++i){ //Fehler3 Wir mussen die Schleife n -mal wiederholen, sizeof(array) ist nur die grosse von Zeiger
                     //nicht die grosse, der reservierten Speicherzellen
       array[i]=1./(i+1); //Fehler4 Wir mussen den Punkt nachdem 1 haben um eine reelle Ergebnis zu bekommen
  }
  printarray_double(array,n);//Fehler5: wir mussen ein neue Funktion schreiben um reelle arrays zu ausdrucken
  free(array);
  return 0;
}
