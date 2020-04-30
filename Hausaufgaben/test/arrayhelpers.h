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
void array_ausgabe_int(int* array, const int length, int option);

/*
*   This function initilizes a list of int data with rising values
*   @length (in): number of datapoints
*
*   @array (out): array with all the data_points.
*/
int* array_init_std_int(int length);

/*
*   This function initilizes an int array with all the same value.
*
*   @length (in): number of datapoints
*   @value (in): value for all datapoints
*/
int* array_init_value_int(int length, int value);

/*
*   This function rotates the array by one
*
*   @array (in): array to be rotated.
*   @length (in): number of datapoints.
*   
*/
void rotate_one_int(int *array, const int length);

/*
*   This function rotates the array by a given amount.
*
*   @array (in): array to be rotated.
*   @length (in): number of datapoints.
*   @k (in): number of rotations to be executed
*/
void rotate_k_int(int *array, const int length, int k);

/*
*   This function turns a given array around.
*
*   @array (in): array to be turned.
*   @length (in): number of datapoints.
*/
void turn_array_int(int *array, const  int length);

int search_array_int(int* array1, int length1, int* array2, int length2);

/*
*   This function swaps two elements of a given array in place.
*
*   @array (in): array to be turned.
*   @length (in): number of datapoints.
*   @k1 (in): position of the first element.
*   @k2 (in): position of the second element.
*/
void swap_pos_int(int* array, int length,  int k1, int k2);

/**
 *  This function creates an array of data points to be used as input values
 *  for numerical simulations or other puroposes. It creates n datapoints
 *  between a and b where a is the first and b the last point.
 * 
 *  @a (in): start of the interval.
 *  @b (in): end of the interval.
 *  @n (in): number of datapoints.
 * 
 *  @data_array (out): array of all data points.
 **/
double* create_array_double(double a, double b, int n);

/**
 *  This function creates data points from an input array and a given function.
 *  The function has to match the data type of the array and the first argument
 *  of the function has to a be the input variable of the function.
 * 
 *  \param data: array of input data points
 *  \param func: pointer of the function used for creating data points
 *  \param length: number of data points
 *  \param p: list of all parameters for the function.
 * 
 *  \return array of data points created by the function.
 * */
double* create_array_function_double(double* data, double(*func)(double, void*), const int length, void *p);

/*
*   This function prints an array of douuble-data to stdout either in line
*   or one data point per line.
*   
*   @array (in): list of data to be printed.
*   @length (in): number of datapoints.
*   @option (in): 1 for line-seperated values, 2 for semicolon-spereated values
*/
void array_ausgabe_double(double* array, const int length, int option);

#endif