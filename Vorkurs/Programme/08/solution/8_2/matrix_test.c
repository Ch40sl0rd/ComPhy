
#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>

void print_usage(void);

int main(int argc, char** argv){
  if( argc < 2 ){
    print_usage();
    exit(MATRIX_INVALID_CMDLINE_ARGUMENT_ERROR);
  }

  unsigned int const n = atoi(argv[1]);

  // allocate some matrices
  matrix_t * a = matrix_alloc(n);
  matrix_t * b = matrix_alloc(n);
  matrix_t * c = matrix_alloc(n);

  // set matrix elements to 
  matrix_set_row_column(a);

  printf("a_ij = j + i*n:\n");
  matrix_print(a);

  // transpose 'a' and store in 'b'
  matrix_transpose(b,a);

  printf("b_ij = a_ji :\n");
  matrix_print(b);

  // perform matrix multiplication a*b and store result in c
  matrix_mult(c,a,b);

  printf("c_ij = a_ik * b_kj:\n");
  matrix_print(c);

  return 0;
}

void print_usage(void){
  printf("Usage: ./matrix_test <matrix size>\n");
}
