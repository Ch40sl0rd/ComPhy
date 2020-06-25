// make
// ./hausaufgabe1
// Lars Döpper
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void print_data2file(char* filename,double *data, int length){
    if(data==NULL){
        printf("compromised data given, abort workflow!\n");
        exit(2);
    }
    FILE *file;
    if(filename==NULL){
        file = stdout;
    }
    else{
        file = fopen(filename, "w");
    }

    if(!file){
        printf("File could not be created!\n");
        exit(3);
    }

    for(int i=0; i<length; i++){
        fprintf(file, "%le\n", data[i]);
    }
    if(file!=stdout){
        fclose(file);
    }
}

void print_data_table_double(char* filename, double *x, double *y, const int length){
    if(!x || !y){
        printf("Can't handle empty pointers.\n");
        exit(404);
    }

    FILE *file;
    if(filename==NULL){
        file = stdout;
    }
    else{
        file = fopen(filename, "w");
    }

    if(!file){
        printf("File could not be created!\n");
        exit(3);
    }

    for(int i=0; i<length; i++){
        fprintf(file, "%15.6e\t%15.6e\n", x[i], y[i]);
    }
    if(file!=stdout){
        fclose(file);
    }
}

void print_table2file(char* const filename, double **table, int rows, int cols){
    FILE *file_out;

    if(filename==NULL){
        file_out = stdout;
    }
    else{
        file_out = fopen(filename, "w");
    }

    if(!file_out){
        printf("File could not be created!\n");
        exit(3);
    }

    for(int i=0; i<rows; i++){
        for(int j=0; j<cols-1; j++){
            fprintf(file_out, "%15.6e\t", table[i][j]);
        }
        fprintf(file_out, "%15.6e\n", table[i][cols-1]);
    }
    if(file_out!=stdout){
        fclose(file_out);
    }
}

void print_data2file_func(char*filename,  double *data, int length , double(*fp)(double)){
    FILE *file_out = fopen(filename, "w");
    if(!file_out){
        printf("Error while openig the file.\n");
        exit(2);
    }
    for(int i=0; i<length; i++){
        fprintf(file_out, "%le\t%le\n", data[i], fp(data[i]));
    }
    fclose(file_out);
}
double* read_data_file(char* const filename, int *length, int binary){
    FILE *file =NULL;
    double* data;
    int reads, error;
    if(filename==NULL){
        printf("Kein Dateiname angegeben.\n");
        return NULL;
    }
    if(binary){
        //reading all the binary data from file
        file = fopen(filename, "rb");
        if(!file){
            printf("Fehler bei der Öffnung der Datei.\n");
            exit(2);
        }
        reads = fread((void*)length, sizeof(int), 1, file);
        if(reads!=1){
            printf("Fehler bei lesen der Länge.\n");
            exit(3);
        }
        printf("Anzahl der Datenpunkte: %d\n", (*length));
        data = (double*)malloc(sizeof(double)*(*length));
        reads = fread( (void*)data, sizeof(double), (*length), file);
        if(reads != (*length)){
            printf("Fehler bei Lesen der eigentlichen Daten.\n");
            exit(3);
        }
        if(ferror(file) !=0){
            printf("ferror of file != 0.\n");
            exit(3);
        }
        error = fclose(file);
        if(error){
            printf("Error while closing the file.\n");
            exit(4);
        }
        return data;

    }
    //sonst normales auslesen der Daten
    file = fopen(filename, "r");
    if(!file){
        printf("Fehler bei der Öffnung der Datei.\n");
        exit(2);
    }

    fscanf(file, "%d\n", length);
    printf("Die Länge der Datei beträgt: %d\n", (*length));
        
    data = (double*)malloc(sizeof(double)*(*length));
    for(int i=0; i<(*length); i++){
        fscanf(file, "%lf\n", &data[i]);
    }
    fclose(file);
    return data;
}