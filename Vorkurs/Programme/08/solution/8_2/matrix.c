#include "matrix.h"
#include "fatal_error.h"

#include <stdlib.h>
#include <stdio.h>
// for memset
#include <string.h>

matrix_t * matrix_alloc(unsigned int const n){
  matrix_t * m = (matrix_t*)malloc(sizeof(matrix_t));
  fatal_error((void*)m == NULL,
              "[matrix_alloc] memory allocation for 'm' failed!\n",
              MATRIX_ALLOC_ERROR);

  m->n = n;

  m->mem = (double*)malloc(n*n*sizeof(double));
  m->M = (double**)malloc(n*sizeof(double*));
  fatal_error((void*)m->mem == NULL || (void*)m->M == NULL,
              "[matrix_alloc] memory allocation for m->mem or m->M failed!\n",
              MATRIX_ALLOC_ERROR);

  for(unsigned int r = 0; r < n; ++r){
    m->M[r] = m->mem + r*n;
  }

  return m;
}

void matrix_free(matrix_t * m){
  fatal_error( (void*)m == NULL,
               "[matrix_free] 'm' cannot be NULL!\n",
               MATRIX_INVALID_ARGUMENT_ERROR );

  free(m->mem);
  free(m->M);
  free(m);
}

matrix_t * matrix_copy(matrix_t const * const m){

  fatal_error((void*)m == NULL,
              "[matrix_copy] 'm' cannot be NULL!\n",
              MATRIX_NULL_POINTER_ERROR);

  matrix_t * m_copy = matrix_alloc(m->n);
  memcpy((void*)m_copy->mem, (void*)m->mem, m->n*m->n*sizeof(double) );
  memcpy((void*)m_copy->M, (void*)m->M, m->n*sizeof(double*) );

  return m_copy;
}

void matrix_transpose(matrix_t * const out, matrix_t const * const in){
  fatal_error((void*)out == NULL ||
              (void*)in == NULL,
              "[matrix_transpose] 'out' or 'in' cannot be NULL!\n",
              MATRIX_NULL_POINTER_ERROR);
  fatal_error(out == in,
              "[matrix_transpose] 'out' and 'in' cannot be the same!\n",
              MATRIX_INVALID_ARGUMENT_ERROR);
  fatal_error(out->n != in->n,
              "[matrix_transpose] 'out' and 'in' must have the same size!\n",
              MATRIX_SIZE_ERROR);

  for(unsigned int r = 0; r < in->n; ++r){
    for(unsigned int c = 0; c < in->n; ++c){
      out->M[r][c] = in->M[c][r];
    }
  }
}

void matrix_transpose_inplace(matrix_t * const m){
  fatal_error( (void*)m == NULL,
               "[matrix_transpose_inplace] 'm' cannot be NULL!\n",
               MATRIX_NULL_POINTER_ERROR );
  
  // we first create a copy of 'm'  
  matrix_t * m_copy = matrix_copy(m);
  
  // now we can perform the transposition
  matrix_transpose(m, m_copy); 

  // now we free the memory that we just used
  matrix_free(m_copy);
}

void matrix_set_id(matrix_t * const m){
  fatal_error( (void*)m == NULL,
               "[matrix_set_id] 'm' cannot be NULL!\n",
               MATRIX_NULL_POINTER_ERROR );

  for(unsigned int r = 0; r < m->n; ++r){
    for(unsigned int c = 0; c < m->n; ++c){
      if(r == c){
        m->M[r][c] = 1.0;
      } else {
        m->M[r][c] = 0.0;
      }
    }
  }
}

void matrix_set_row_column(matrix_t * const m){
  fatal_error( (void*)m == NULL,
               "[matrix_set_row_column] 'm' cannot be NULL!\n",
               MATRIX_NULL_POINTER_ERROR );

  for(unsigned int r = 0; r < m->n; ++r){
    for(unsigned int c = 0; c < m->n; ++c){
      m->M[r][c] = c + m->n*r;
    }
  }

}

void matrix_print(matrix_t const * const m){
  fatal_error( (void*)m == NULL,
               "[matrix_print] 'm' cannot be NULL!\n",
               MATRIX_NULL_POINTER_ERROR );

  printf("\n");
  for(unsigned int r = 0; r < m->n; ++r){
    for(unsigned int c = 0; c < m->n; ++c){
      printf("%8.2f  ", m->M[r][c]);
    }
    printf("\n");
  }
  printf("\n");
}

void matrix_mult(matrix_t * const c, matrix_t const * const a, matrix_t const * const b){
  fatal_error((void*)a == NULL ||
              (void*)b == NULL ||
              (void*)c == NULL,
              "[matrix_mult] Neither 'a' nor 'b' nor 'c' can be NULL!\n",
              MATRIX_NULL_POINTER_ERROR );
  fatal_error(a == c || b == c,
              "[matrix_mult] Neither 'a' nor 'b' can be the same as 'c'!\n",
              MATRIX_INVALID_ARGUMENT_ERROR);
  fatal_error(a->n != b->n ||
              b->n != c->n ||
              a->n != c->n,
              "[matrix_mult] Matrices 'a', 'b' and 'c' must have the same size!\n",
              MATRIX_SIZE_ERROR);


  // we first zero out 'c'
  memset((void*)c->mem, 0, sizeof(double)*c->n*c->n);

  // row index for 'c'
  for(unsigned int i = 0; i < c->n; ++i){
    // column index for 'c'
    for(unsigned int j = 0; j < c->n; ++j){
      // row index for 'b'
      for(unsigned int k = 0; k < c->n; ++k){
        // C_ij = \sum_k A_ik * B_kj
        c->M[i][j] += a->M[i][k] * b->M[k][j];
      }
    }
  }
}
