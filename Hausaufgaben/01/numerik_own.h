// make
// ./hausaufgabe1
// Lars Döpper
#ifndef NUMERIK_OWN_H
#define NUMERIK_OWN_H

/**
*   This function calculates the first symmetric derivation of a given function.
*   \param x (in): point to calculate the value at
*   \param h (in): distance for the derivation
*   \param func (in): function to calculate the derivation to.
*
*   \return: value of the first derivation at position x.
*/
double derivate_sym_one(double x, double h, double(*func)(double));

/**
 *  This function calculates the first symmetric derivation from a
 *  array of data points (y-coord.) and returns an array of 
 *  y-coord for the first derivation.
 * 
 *  \param x: array of data points for input purposes
 *  \param length: number of data points
 *  \param a: x value of the first data point
 *  \param b: x value of the last data point.
 * 
 *  \return: array of data points for the fist derivation
 * */
double* derivate_sym_one_array(double* data, double a, double b, const int length);

/**
*   This function calculates the second symmetric derivation of a given function.
*
*   \param x (in): point to calculate the value at
*   \param h (in): distance for the derivation
*   \param func (in): function to calculate the derivation to.
*
*   \return: value of the second derivation at position x.
**/
double derivate_sym_two(double x, double h, double(*func)(double));

/**
 *  This function calculates the second symmetric derivation from a
 *  array of data points (y-coord.) and returns an array of 
 *  y-coord for the first derivation.
 * 
 *  \param x: array of data points for input purposes
 *  \param length: number of data points
 *  \param a: x value of the first data point
 *  \param b: x value of the last data point.
 * 
 *  \return: array of data points for the fist derivation
 * */
double* derivate_sym_two_array(double* data, double a, double b, const int length);

/**
*   This function calculates the third symmetric derivation of a given function.
*   \param x (in): point to calculate the value at
*   \param h (in): distance for the derivation
*   \param func (in): function to calculate the derivation to.
*
*   \return: value of the third derivation at position x.
*/
double derivate_sym_three(double x, double h, double(*func)(double));

/**
*   This function calculates the position of the zero crossings of
*   a given function by using newton's methode with unknown derivations.
*   
*   \param func : function to use.
*   \param a : first starting position.
*   \param b : second starting position.
*   \param acc : accuracy of the zero corssing position.
*
*   \return: position of the zero crossing within accuracy
*/
double zero_crossing(double (*func)(double, void*), const double x0, const double x1, double acc, void *p);

/**
*   This function will calculate the finite integrale of the given function
*   by using the trapez methode with a defined stepsize.
*
*   \param f (in): function to be integrated.
*   \param a (in): start of the integrale.
*   \param b (in): end of the integrale.
*   \param h (in): stepsize.
*
*   \returns: value of the integrale.
*/
double integrate_trapez(double (*f)(double), double a, double b, double h);

/**
*   This function will calculate the finite integrale of the given function
*   by using the trapez methode with an adaptive stepsize.
*
*   \param f (in): function to be integrated.
*   \param a (in): start of the integrale.
*   \param b (in): end of the integrale.
*   \param acc (in): desired accuracy of the value.
*
*   \returns: value of the integrale.
*/
double integrate_trapez_adap( double(*f)(double), double a, double b, double acc);

/**
 *  This function uses the romberg integration methode to calculate an finite
 *  integral from a to b of the function func. It repeats the process until the
 *  desired accucary is reached.
 * 
 *  \param a: start of the integration interval.
 *  \param b: end of the integration interval.
 *  \param n_0: number of steps in the first iteration.
 *  \param func: function to be integrated.
 *  \param acc: desired accuracy of the integration.
 *  \param m_max: maximum number of steps.
 *  \param p: list of all parameters for the function.
 * 
 *  \return: valuue of the integrale
 * */
double integrate_romberg(double a, double b, int n_0, double (*func)(double, void*), double acc, int m_max, void* p);

/**
    This function represents the integrand for the gamma function

    \param t (in): point t of the gamma function
    \param p (in): void pointer where the first element has to z for gamma(z)

    \return: value of the gamma integrant at t.
*/
double  gamma_integrand(double t, void *p);

/**
 *  This function calculates the value of the gamma function at a point z
 *  where z has to be a real value. 
 * 
 *  \param z (in): point for the gamma function.
 * 
 *  \return: value of the gamma function.
 */
double gamma_func(double z);

/**
 *  This function calculates the points and weights for the gaus-legendre integration.
 * 
 *  \param z (in): point for the gamma function.
 * 
 *  \return: value of the gamma function.
 */
void gaus_legendre(double a, double b, double *x_vals, double *w_vals, int n);

double gl_integrate(double a, double b, double *x, double *w, int n, double(*f)(double, void *), void *p);
#endif