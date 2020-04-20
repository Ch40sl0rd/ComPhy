#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include"arrayhelpers.h"

double abs_d(double x){
    if(x<0){
        return -x;
    }
    return x;
}

int main(){
    //Festlegen  des Dateinamens
    char const * dateiname="cos.txt";
    //Öffnen der Datei im Modus read
    FILE* file = fopen(dateiname, "r");
    //Fehlerüberprüfung
    if(file  == NULL){
        printf("Error bei Dateiöffnen!\n");
        return 1;
    }
    //Auslesen  der Nummer an Zeilen.
    int line_number;
    float x,  cos_x, diff;
    //Liest die Erste Zahl in der Ersten Zeile des Dokuments aus.
    fscanf(file, "%d\n", &line_number);
    //Gibt die Anzahl an Zeilen aus.
    printf("%d\n", line_number);
    //Überprüfung aller nachfolgenden Zeilen der Datei.
    for(int i=1; i<=line_number; i++){
        fscanf(file, "cos(%f) = %f\n", &x, &cos_x);
        diff = cos_x - cos(x);
        if(abs_d(diff)>1e-3){
            printf("In Zeile %d beträgt die Differenz %f\n", i, diff);
        }
    }
    //Schließt die Datei wieder.
    fclose(file);
    return 0;
}