// Copyright © 2017 Martin Ueding <dev@martin-ueding.de>
// Licensed under the MIT license.

// Implementation of the “Game of Life”.
//
// This uses a linearized 2D array with an index function (`idx`).

#include <stdio.h>
#include <stdlib.h>

// For `memset`.
#include <string.h>

// For `usleep`. If usleep cannot be found, switch to sleep instad
// below.
#include <unistd.h>

/**
  Computes the linarized index with periodic boundary conditions.

  One can insert negative indices like `-1` and it will be wrapped around
  periodically. Also one can insert `i = i_max` and it will wrap it around to
  `i = 0`. This way one does not need to bother with the boundary conditions in
  the actual evolution code in `evolve_field`.

  \param[in] i Site along the first dimension.
  \param[in] j Site along the first dimension.
  \param[in] i_max Extent of the array in the first dimension.
  \param[in] j_max Extent of the array in the second dimension.

  \returns Linearized index between `0` and `i_max * j_max - 1`.
  */
int
idx(int const i, int const j, int const i_max, int const j_max) {
    // Apply periodic boundary conditions.
    int const i_wrap = (i + i_max) % i_max;
    int const j_wrap = (j + j_max) % j_max;

    return i_wrap * j_max + j_wrap;
}

/**
  Swaps the contents of two pointers.
  */
void swap(char **a, char **b) {
    char *temp;
    temp = *a;
    *a = *b;
    *b = temp;
}

/**
  Creates a new game field and sets every site to “dead”.

  \param[in] i_max Extent of the array in the first dimension.
  \param[in] j_max Extent of the array in the second dimension.
  */
char *make_field(int const i_max, int const j_max) {
    // Calculate the needed number of bytes for the playing field.
    int const bytes = i_max * j_max * sizeof(char);

    // Allocate the memory and check for errors.
    char *field = malloc(bytes);
    if (field == NULL) {
        abort();
    }

    // Fill the field with zeros, “dead” cells.
    memset(field, 0, bytes);

    return field;
}

/**
  Evolve the game field one time step.

  It is important to have two fields. Otherwise one would start to overwrite the
  one field with new value. One cannot properly count the number of neighbors
  any more. Therefore this function uses two fields and does not overwrite any
  data in `field_in`.

  \param[in] field_in Current state of the game field.
  \param[out] field_out Next state of the game field.
  \param[in] i_max Extent of the array in the first dimension.
  \param[in] j_max Extent of the array in the second dimension.
  */
void evolve_field(char *field_in,
                  char *field_out,
                  int const i_max,
                  int const j_max) {
    // Iterate through the field vertically ...
    for (int i = 0; i < i_max; ++i) {
        // ... and horizontally.
        for (int j = 0; j < j_max; ++j) {
            // Read the old value of the cell.
            char const old_value = field_in[idx(i, j, i_max, j_max)];

            // Count the neighbors. This uses that a dead cell has value `0`, a
            // living cell has value `1`. Therefore the sum of the fields is
            // just the number of the neighbors.
            char const neighbors = field_in[idx(i + 1, j + 1, i_max, j_max)] +
                                   field_in[idx(i + 1, j + 0, i_max, j_max)] +
                                   field_in[idx(i + 1, j - 1, i_max, j_max)] +
                                   field_in[idx(i + 0, j + 1, i_max, j_max)] +
                                   field_in[idx(i + 0, j - 1, i_max, j_max)] +
                                   field_in[idx(i - 1, j + 1, i_max, j_max)] +
                                   field_in[idx(i - 1, j + 0, i_max, j_max)] +
                                   field_in[idx(i - 1, j - 1, i_max, j_max)];

            // Compute the new value according to the rules of the game.
            char new_value;
            if (old_value == 1) {
                if (neighbors < 2) {
                    new_value = 0;
                } else if (neighbors <= 3) {
                    new_value = 1;
                } else {
                    new_value = 0;
                }
            } else {
                if (neighbors == 3) {
                    new_value = 1;
                } else {
                    new_value = 0;
                }
            }

            // Set the new value into the output field.
            field_out[idx(i, j, i_max, j_max)] = new_value;
        }
    }
}

/**
  Prints out the field to the console.

  \param[in] field Game field to print.
  \param[in] i_max Extent of the array in the first dimension.
  \param[in] j_max Extent of the array in the second dimension.
  */
void print_field(char *field, int const i_max, int const j_max) {
    printf("\n\n");
    for (int i = 0; i < i_max; ++i) {
        for (int j = 0; j < j_max; ++j) {
            char value = field[idx(i, j, i_max, j_max)];
            if (value == 0) {
                printf(".");
            } else {
                printf("#");
            }
        }
        printf("\n");
    }
}

int main(int argc, char **argv) {
    int const i_max = 10;
    int const j_max = 20;

    // Create one input and one output field.
    char *field1 = make_field(i_max, j_max);
    char *field2 = make_field(i_max, j_max);

    // Insert the “glider” structure.
    field1[idx(2, 1, i_max, j_max)] = 1;
    field1[idx(3, 2, i_max, j_max)] = 1;
    field1[idx(1, 3, i_max, j_max)] = 1;
    field1[idx(2, 3, i_max, j_max)] = 1;
    field1[idx(3, 3, i_max, j_max)] = 1;

    for (int iteration = 0; iteration < 100; ++iteration) {
        print_field(field1, i_max, j_max);

        // Evolve `field` to `field2`.
        evolve_field(field1, field2, i_max, j_max);

        // Swap the two pointers (not the data!) such that `field1` points to
        // the new state.
        swap(&field1, &field2);

        // Sleep a bit (500 000 micro seconds) before doing the next iteration.
        usleep(500000);
        // if usleep cannot be found, use sleep instead
        // sleep(1);
    }

    // Clean up.
    free(field1);
    free(field2);
}
