//Lars Döpper
//make
// ./hausaufgabe4
#ifndef FILEHELPERS_h
#define FILEHELPERS_H


/**
*   This function prints a dataset to a file or stdout if filename ist NULL
*
*   \param filename (in): name of the file to write to or NULL for stdout
*   \param data (in): data to be printed to the destination.
*   \param length (in): number of datapoints to be printed.
*/
void print_data2file(char* filename,double *data, int length);

/**
*   This function prints pairs of data to a given outpuit-file. The data pairs
*   consist of (x, fp(x))
*
*   \param filename (in): name of the output file.
*   \param data (in): data to be printed.
*   \param length (in): number of data points.
*   \param fp (in): function to be applied to the data.
*/
void print_data2file_func(char*filename,  double *data, int length , double(*fp)(double));

/**
    This function creates a data array form a given file list, the length
    of the data to be importet must be known

    \param filename (in): name of the file containing the data.
    \param length (in): pointer to length or number of data-points in the file.

    \return list (out): list of all data points.
*/
double* read_data_file(char* const filename, int *length, int binary);

/**
 *  This function prints two arrays as a data table with x\t y
 *  for easy data transfer to graphical tools or other programs.
 * 
 *  \param filename: name of the file to be printed to.
 *  \param x: array of x values:
 *  \param y: array of y values.
 *  \param length: nuumber of datapoints.
 * 
 * */
void print_data_table_double(char* filename, double *x, double *y, const int length);

void print_table2file(char* const filename, double **table, int rows, int cols);
#endif