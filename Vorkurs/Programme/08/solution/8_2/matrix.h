#ifndef MATRIX_H
#define MATRIX_H

/**
 * @brief enumeration for error conditions
 *
 * An enumeration is a special kind of mapping between names and integer
 * values which can be used much like pre-processor constants. The added
 * benefit is that the data type thus defined, matrix_error_t, is an actual
 * data type and variables of this type can only take the values below
 * and no other.
 */
typedef enum matrix_error_t { 
  MATRIX_ALLOC_ERROR = 128,
  MATRIX_SIZE_ERROR = 129,
  MATRIX_INVALID_ARGUMENT_ERROR = 130,
  MATRIX_NULL_POINTER_ERROR = 131,
  MATRIX_INVALID_CMDLINE_ARGUMENT_ERROR = 132
} matrix_error_t;

/**
 * @brief data type to store an n*n matrix
 */
typedef struct matrix_t {
  double **M;
  double *mem;
  unsigned int n;
} matrix_t;

/**
 * @brief Allocator for n*n (double) matrix
 *
 * Upon failure, exit will be called.
 *
 * @param n Matrix size: n^2
 *
 * @return Pointer to newly allocated matrix. The user is responsible for
 *         deallocation via 'matrix_free'.
 */
matrix_t * matrix_alloc(unsigned int const n);

/**
 * @brief Create a copy of an n*n (double) matrix
 *
 * @param m Pointer to matrix allocated with 'matrix_alloc'. 'm' is
 *          not modified!
 *
 * @return Pointer to a newly allocated matrix, all elements of which
 *         have been set equal to the original matrix 'm'. 
 */
matrix_t * matrix_copy(matrix_t const * const m);

/**
 * @brief Free memory allocated via 'matrix_alloc'
 *
 * @param m Pointer to a matrix allocated with 'matrix_alloc'.
 */
void matrix_free(matrix_t * m);

/**
 * @brief Transpose a matrix "in place"
 *
 * A copy of the matrix 'm' is made and then transposed.
 *
 * @param m Matrix to be transposed as well as the output matrix.
 */
void matrix_transpose_inplace(matrix_t * const m);


/**
 * @brief Transpose a matrix
 *
 * 'out' and 'in' must have the same size and cannot be the
 * same objects
 *
 * @param out A matrix allocated via 'matrix_alloc' that the result
 *            of the transposition should be stored in. 
 * @param in A matrix allocated via 'matrix_alloc' which should
 *           be transposed.
 */
void matrix_transpose(matrix_t * const out, matrix_t const * const in);

/**
 * @brief Set a matrix to the identity
 *
 * The diagonal elements of 'm' are set to 1.0 while all other
 * elements are set to 0.0
 *
 * @param m Matrix allocated with 'matrix_alloc'.
 */
void matrix_set_id(matrix_t * const m);

/**
 * @brief Set matrix elements to their linear index
 *
 *  Each matrix element m_ij will be set to j + n*i, cast to double.
 *
 * @param m Matrix allocated with 'matrix_alloc'.
 */
void matrix_set_row_column(matrix_t * const m);

/**
 * @brief Perform the matrix multiplication c = a*b
 *
 * Multiplies matrices a and b and stores the result in c.
 * a,b and c must have the same dimensions
 * a and b can be the same, but neither a nor b can be the same as c
 *
 *
 * @param c Matrix allocated via 'matrix_alloc' for output
 * @param a left matrix to be post-multiplied by 'b'
 * @param b right matrix to be pre-multiplied by 'a'
 */
void matrix_mult(matrix_t * const c, matrix_t const * const a, matrix_t const * const b);

/**
 * @brief Print matrix elements to screen
 *
 * @param m Matrix allocated via 'matrix_alloc' for output
 */
void matrix_print(matrix_t const * const m);

#endif // MATRIX_H
