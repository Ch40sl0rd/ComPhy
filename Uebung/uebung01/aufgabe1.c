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

int main(){
    char* filename1 = "dataset1";
    char *filename2 = "dataset2";
    char* file_ausgabe = "ausgabe.txt";
    int length1, length2;

    double *list_data = read_data_file(filename1, &length1, 0);
    double *list_bin_data = read_data_file(filename2, &length2, 1);
    double *list_sum = sum_sqr_sum(list_data, length1);

    print_data_file(NULL, list_data, length1);
    print_data_file(NULL, list_bin_data, length2);

    print_data_file(file_ausgabe, list_sum, length1);

    print_data2file_func("ausgabe2.txt", list_bin_data, length2, sin);

    free(list_data);
    free(list_bin_data);
}