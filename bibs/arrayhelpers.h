#ifndef ARRAYHELPERS_H
#define ARRAYHELPERS_H

/*
*   This function prints an array of int-data to stdout either in line
*   or one data point per line.
*   
*   @array (in): list of data to be printed.
*   @length (in): number of datapoints.
*   @option (in): 1 for line-seperated values, 2 for semicolon-spereated values
*/
void array_ausgabe(int* array, const int length, int option);

/*
*   This function initilizes a list of int data with rising values
*   @length (in): number of datapoints
*
*   @array (out): array with all the data_points.
*/
int* array_init_std(int length);

/*
*   This function initilizes an int array with all the same value.
*
*   @length (in): number of datapoints
*   @value (in): value for all datapoints
*/
int* array_init_value(int length, int value);

/*
*   This function rotates the array by one
*
*   @array (in): array to be rotated.
*   @length (in): number of datapoints.
*   
*/
void rotate_one(int *array, const int length);

/*
*   This function rotates the array by a given amount.
*
*   @array (in): array to be rotated.
*   @length (in): number of datapoints.
*   @k (in): number of rotations to be executed
*/
void rotate_k(int *array, const int length, int k);

/*
*   This function turns a given array around.
*
*   @array (in): array to be turned.
*   @length (in): number of datapoints.
*/
void turn_array(int *array, const  int length);

int search_array(int* array1, int length1, int* array2, int length2);

/*
*   This function swaps two elements of a given array in place.
*
*   @array (in): array to be turned.
*   @length (in): number of datapoints.
*   @k1 (in): position of the first element.
*   @k2 (in): position of the second element.
*/
void swap_pos(int* array, int length,  int k1, int k2);

#endif