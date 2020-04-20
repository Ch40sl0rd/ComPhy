#include <stdio.h>

void array_ausgabe(int array[], const int length, int option){
    if(option == 1){
        for(int i = 0; i < length; i++){
            printf("%d\n", array[i]);
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

void array_init_std(int* array, int length){
    for(int i = 0; i<length; i++){
        array[i] = i+1;
    }
}

void array_init_value(int *array, int length, int value){
    for(int i=0; i<length; i++){
        array[i] = value;
    }
}

void rotate_1(int *array, const int length){
    int temp = array[0];
    for(int i = 0;  i< length-1; i++){
        array[i] = array[i+1];
    }
    array[length-1] = temp;
}

void rotate_k(int *array, const int length, int k){
    for(int i=0; i<k; i++){
        rotate_1(array, length);
    }
}

void turn_array(int *array, const  int length){
    int  temp;
    for (int i=0; i<length/2.;i++){
        temp = array[i];
        array[i] = array[length-1-i];
        array[length-1-i] =  temp;
    }
}

int search_array(int* array1, int length1, int* array2, int length2){
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

void swap_pos(int* array, int length,  int k1, int k2){
    int temp = array[k1];
    array[k1]  =  array[k2];
    array[k2] = temp;
}