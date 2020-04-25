#ifndef FATAL_ERROR_H
#define FATAL_ERROR_H

/**
 * @brief Function to terminate program with signal
 *
 * If 'condition' evaluates non-zero (true), fatal_error will print
 * the string 'msg' to stdout and terminate the program using exit,
 * sending the signal passed via 'signal'.
 *
 * @param condition Boolean condition which causes program termination
 *                  if it evaluates non-zero.
 * @param msg Message to be printed to stdout, should be a null-terminated
 *            string.
 * @param signal Signal to be passed to exit.
 */
void fatal_error(int const condition, char const * const msg, int const signal); 

#endif // FATAL_ERROR_H
