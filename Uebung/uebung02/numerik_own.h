#ifndef NUMERIK_OWN_H
#define NUMERIK_OWN_H

/*
*   This function calculates the first symmetric derivation of a given function.
*   x (in): point to calculate the value at
*   h (in): distance for the derivation
*   func (in): function to calculate the derivation to.
*
*   return: value of the first derivation at position x.
*/
double derivate_sym_one(double x, double h, double(*func)(double, void*), void* p);

/*
*   This function calculates the second symmetric derivation of a given function.
*   x (in): point to calculate the value at
*   h (in): distance for the derivation
*   func (in): function to calculate the derivation to.
*
*   return: value of the second derivation at position x.
*/
double derivate_sym_two(double x, double h, double(*func)(double, void*), void* p);

/*
*   This function calculates the third symmetric derivation of a given function.
*   x (in): point to calculate the value at
*   h (in): distance for the derivation
*   func (in): function to calculate the derivation to.
*
*   return: value of the third derivation at position x.
*/
double derivate_sym_three(double x, double h, double(*func)(double, void*), void* p);

/*
*   This function calculates the position of the zero crossings of
*   a given function by using newton's methode with unknown derivations.
*   
*   @func (in): function to use.
*   @a (in): first starting position.
*   @b (in): second starting position.
*   @acc (in): accuracy of the zero corssing position.
*
*   @return: position of the zero crossing within accuracy
*/
double zero_crossing_newton(double (*func)(double), double x0, double x1, double acc);

/*
*   This function will calculate the finite integrale of the given function
*   by using the trapez methode with a defined stepsize.
*
*   @f (in): function to be integrated.
*   @a (in): start of the integrale.
*   @b (in): end of the integrale.
*   @h (in): stepsize.
*
*   @returns: value of the integrale.
*/
double integrate_trapez(double (*f)(double), double a, double b, double h);

/*
*   This function will calculate the finite integrale of the given function
*   by using the trapez methode with an adaptive stepsize.
*
*   @f (in): function to be integrated.
*   @a (in): start of the integrale.
*   @b (in): end of the integrale.
*   @acc (in): desired accuracy of the value.
*
*   @returns: value of the integrale.
*/
double integrate_trapez_adap( double(*f)(double), double a, double b, double acc);

/*
    This function represents the integrand for the gamma function

    @t (in): point t of the gamma function
    @p (in): void pointer where the first element has to z for gamma(z)

    @return: value of the gamma integrant at t.
*/
double  gamma_integrand(double t, void *p);

/*
 *  This function calculates the value of the gamma function at a point z
 *  where z has to be a real value. 
 * 
 *  @z (in): point for the gamma function.
 * 
 *  @return: value of the gamma function.
 */
double gamma_func(double z);
#endif