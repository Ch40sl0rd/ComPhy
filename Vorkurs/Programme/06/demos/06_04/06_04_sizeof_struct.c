#include <stdio.h>

typedef struct test_t {
  char zeichen;
  int zahl;
  double gleitzahl;
} test_t;

int main(void){
  printf("sizeof(test_t) = %lu\n"
         "sizeof(char) = %lu + sizeof(int) = %lu + sizeof(double) = %lu\n",
         sizeof(test_t), 
         sizeof(char), sizeof(int), sizeof(double));
  return 0;
}
