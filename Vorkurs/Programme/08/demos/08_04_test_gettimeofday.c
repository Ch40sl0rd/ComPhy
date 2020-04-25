#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>
#include <stdlib.h>
int main(void){
  struct timeval anfang, ende;
  gettimeofday(&anfang, NULL);
  sleep(1);
  gettimeofday(&ende, NULL);
  
  double sec = (double)(ende.tv_sec-anfang.tv_sec) + 1.0e-6*(double)(ende.tv_usec-anfang.tv_usec);
  printf("Es sind %lf Sekunden vergangen\n", sec);
  
  return 0;
} 

  
