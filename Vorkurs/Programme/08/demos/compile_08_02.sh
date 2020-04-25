gcc -c -I/usr/include/gsl/ 08_02_test_gsl_rng.c
gcc -o 08_02_test_gsl_rng 08_02_test_gsl_rng.o -lgsl -lgslcblas -L/usr/lib/x86_64-linux-gnu/
