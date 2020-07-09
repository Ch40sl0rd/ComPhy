#ifndef RANDOM_H
#define RANDOM_H
//function to initialize random number generator
void init_genrand(unsigned long s);

//function for generating a random number between 0 and 1
double genrand_res53(void);
#endif