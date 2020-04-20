#include <stdio.h>
//Fahlende include von stdlib
#include <stdlib.h>
#include "arrayhelpers.h"

int main(){
    double *array;
    int n,i;
    
    n=42;
    //Falsche Größe für malloc, benötigt sizeof Datentyp
    array = (double*)malloc(n*sizeof(double));
    if(array==NULL){
        printf("Fehler bei Speicherallkokierung\n");
        //Programm wird ohne Return nicht abgebrochen.
        return -1;
    }
    for(i=0; i<n; i++){
        //Fehler mit Integer Division
        array[i] = 1./(i+1);
    }
    //Arrayausgabe funktioniert eigentlich nur mit int Werten
    array_ausgabe(array, n, 2);
    free(array);
    return 0;
}