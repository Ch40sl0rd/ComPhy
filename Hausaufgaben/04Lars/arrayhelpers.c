//Lars Döpper
//make
// ./hausaufgabe4
#include <stdio.h>
#include <stdlib.h>

void array_ausgabe_int(int* array, const int length, int option){
    if(option == 1){
        for(int i = 0; i < length; i++){
            printf("%d \n", array[i]);
        }
    }
    else if(option == 2){
        for(int i = 0; i < length-1; i++){
            printf("%d;", array[i]);
        }
        printf("%d\n", array[length-1]);
    }
    else{
        printf("ERROR. Bitte wähle eine mögliche Ausgabe mit  1 für Zeilenweise und 2 für komma-getrennte Werte.\n");
    }
}

int* array_init_std_int(int length){
    int* data = (int*)malloc(sizeof(int)*length);
    if(!data){
        printf("Error while allocating storage.\n");
        exit(1);
    }
    for(int i=0; i<length; i++){
        data[i] = i+1;
    }
    return data;
}

int* array_init_value_int(int length, int value){
    int* data = array_init_std_int(length);
    for(int i=0; i<length; i++){
        data[i] = value;
    }
    return data;
}

void rotate_one_int(int *array, const int length){
    int temp = array[0];
    for(int i = 0;  i< length-1; i++){
        array[i] = array[i+1];
    }
    array[length-1] = temp;
}

void rotate_k_int(int *array, const int length, int k){
    for(int i=0; i<k; i++){
        rotate_one_int(array, length);
    }
}

void turn_array_int(int *array, const  int length){
    int  temp;
    for (int i=0; i<length/2.;i++){
        temp = array[i];
        array[i] = array[length-1-i];
        array[length-1-i] =  temp;
    }
}

int search_array_int(int* array1, int length1, int* array2, int length2){
    //This function searches for the position of array 1 in array 2.
    int n = 0;
    while(array2[0] != array1[n]){
        n++;
        if( (n-1) == length1){
            return -1;
        }
    }
    for(int i = 0; i<length2; i++){
        if(array2[i] != array1[i+n]){
            return -1;
        }
    }

    return n;
}

void swap_pos_int(int* array, int length,  int k1, int k2){
    int temp = array[k1];
    array[k1]  =  array[k2];
    array[k2] = temp;
}

double* create_array_double(double a, double b, int n){
    double* data_ret;
    double step;

    data_ret = (double*)malloc(sizeof(double)*n);
    if(!data_ret){
        printf("Error while assigning local storage. Programm will be terminated.\n");
        exit(1);
    }

    step = (b-a)/(double)(n-1);
    for(int i=0; i<n; i++){
        data_ret[i] = a+ i*step;
    }
    return data_ret;
}

double* create_array_function_double(double* data, double(*func)(double, void*), const int length, void *p){
    if(data == NULL){
        printf("Error while handling an empty data array.\n");
        exit(404);
    }

    double * data_ret;
    data_ret = (double*)malloc(sizeof(double)*length);
    if(!data_ret){
        printf("Error while allocating local memory.\n");
        exit(1);
    }
    for(int i=0; i<length; i++){
        data_ret[i] = func(data[i], p);
    }
    return data_ret;
}

void array_ausgabe_double(double* array, const int length, int option){
    const char* format;
    if(array == NULL){
        printf("Error while handling an empty array.\n");
        exit(404);
    }
    if(option == 1){
        format = "%15.6e\n";
    }
    else if(option == 2){
        format = "%15.6e ; ";
    }
    for(int i=0; i<length; i++){
        printf(format, array[i]);
    }
    printf("\n");
}

double** create_2d_array(int rows, int cols){
    double *values, **ret_array;
    //allocate 2d array to save values for each step
    values = (double*)malloc(sizeof(double)*cols*rows);
    ret_array = (double**)malloc(sizeof(double*)*rows);
    if(!values || !ret_array){
        printf("[create_2d_array Error while allocating memory.\n");
        abort();
    }
    for(int i=0; i<rows; i++){
        ret_array[i] = values+(i*cols);
    }
    return ret_array;
}

void free2d(double **array){
    if(array==NULL){
        printf("[free_2d] Can#t free NULL-array.\n");
        return;
    }
    free(array[0]); free(array);
}