#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
    This function creates a data array form a given file list, the length
    of the data to be importet must be known

    @filename (in): name of the file containing the data.
    @length (in): pointer to length or number of data-points in the file.

    @list (out): list of all data points.
*/
double* read_data_file(char* const filename, int *length, int binary){
    if(filename==NULL){
        printf("Kein Dateiname angegeben.\n");
        return NULL;
    }
    FILE *file = fopen(filename, "r");
    if(!file){
        printf("Fehler bei der Öffnung der Datei.\n");
        exit(2);
    }

    fscanf(file, "%d\n", length);
    printf("Die Länge der Datei beträgt: %d\n", (*length));
        
    double *list = (double*)malloc(sizeof(double)*(*length));
    for(int i=0; i<(*length); i++){
        fscanf(file, "%lf\n", &list[i]);
    }
    fclose(file);
    return list;
}

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
*   This function creates a sum of the squares of the given data and returns
*   a list of all values.
*   @data (in): array of data of type double to be worked with.
*   @length (in): length of the data list.
*
*   @ret (out): list of partial summs.
*/
double* sum_sqr_sum(double* array, int length){
    double *ret = (double*)malloc(sizeof(double)*length);
    if(!ret){
        printf("Fehler bei Speicherallokierung!\n");
        exit(1);
    }
    double sum=0.;
    for(int i=0; i<length; i++){
        sum += array[i]*array[i];
        ret[i] = sum;
    }
    return ret;
}

int main(){
    char* filename = "dataset1";
    char* file_ausgabe = "ausgabe.txt";
    FILE *file_in = fopen(filename, "r");
    int length;

    double *list_data = read_data_file(filename, &length, 0);
    print_data_file(NULL, list_data, length);
    double *list_sum = sum_sqr_sum(list_data, length);
    print_data_file(file_ausgabe, list_sum, length);

    free(list_data);
    fclose(file_in);
}