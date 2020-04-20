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
    float x, cos_x, diff;
    //Öffnen der Datei im Modus read
    FILE* file = fopen(dateiname, "r");
    //Fehlerüberprüfung
    if(file  == NULL){
        printf("Error bei Dateiöffnen!\n");
        return 1;
    }
    //Gehe die Datei Zeile für Zeile durch und merke dir die Zeilennummer
    int i=0;
    //feof() gibt zurück, ob der Dateipointer schon am Ende angekommen ist.
    while(!feof(file)){
        fscanf(file, "cos(%f) = %f\n", &x, &cos_x);
        //printf("Test  %d\n", i);
        i++;
    }
    int line_number = i;
    //Gibt die Anzahl an Zeilen aus.
    printf("Es gibt %d Zeilen\n", line_number);
    //Pointer an den Anfang zurücksetzen
    rewind(file);
    //Allokierung von Speicher für Wertepaare
    float *x_a, *cos_x_a, *cos_dat_x_a;
    x_a = (float*)malloc(sizeof(float)*line_number);
    cos_x_a = (float*)malloc(sizeof(float)*line_number);
    cos_dat_x_a  = (float*)malloc(sizeof(float)*line_number);
    //Überprüfen der Speicherallokierung
    if( !x_a || !cos_x_a || !cos_dat_x_a){
        printf("Error bei Speicherzuweisung!\n");
        exit(1);
    }



    //Überprüfung aller Zeilen der Datei und speicherung der Werte.
    for(i=0; i<line_number; i++){
        fscanf(file, "cos(%f) = %f\n", &x, &cos_x);
        diff = cos_x - cos(x);
        if(abs_d(diff)>1e-3){
            printf("In Zeile %d beträgt die Differenz %f\n", i, diff);
        }
        x_a[i] = x;
        cos_x_a[i] = (float)cos(x);
        cos_dat_x_a[i] = cos_x;
    }
    //Schließt die Datei wieder.
    free(x_a);
    free(cos_x_a);
    free(cos_dat_x_a);
    fclose(file);
    return 0;
}