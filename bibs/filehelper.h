#include <stdio.h>
#include <stdlib.h>

/*
*   This function prints a dataset to a file or stdout if filename ist NULL
*
*   @filename (in): name of the file to write to or NULL for stdout
*   @data (in): data to be printed to the destination.
*   @length (in): number of datapoints to be printed.
*/
void print_data_file(char* filename,double *data, int length){
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

/*
*   This function prints pairs of data to a given outpuit-file. The data pairs
*   consist of (x, fp(x))
*
*   @filename (in): name of the output file.
*   @data (in): data to be printed.
*   @length (in): number of data points.
*   @fp (in): function to be applied to the data.
*/
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

/*
    This function creates a data array form a given file list, the length
    of the data to be importet must be known

    @filename (in): name of the file containing the data.
    @length (in): pointer to length or number of data-points in the file.

    @list (out): list of all data points.
*/
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