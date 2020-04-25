#ifndef FILEHELPERS_h
#define FILEHELPERS_H


/*
*   This function prints a dataset to a file or stdout if filename ist NULL
*
*   @filename (in): name of the file to write to or NULL for stdout
*   @data (in): data to be printed to the destination.
*   @length (in): number of datapoints to be printed.
*/
void print_data_file(char* filename,double *data, int length);

/*
*   This function prints pairs of data to a given outpuit-file. The data pairs
*   consist of (x, fp(x))
*
*   @filename (in): name of the output file.
*   @data (in): data to be printed.
*   @length (in): number of data points.
*   @fp (in): function to be applied to the data.
*/
void print_data2file_func(char*filename,  double *data, int length , double(*fp)(double));

/*
    This function creates a data array form a given file list, the length
    of the data to be importet must be known

    @filename (in): name of the file containing the data.
    @length (in): pointer to length or number of data-points in the file.

    @list (out): list of all data points.
*/
double* read_data_file(char* const filename, int *length, int binary);
#endif