#include "fatal_error.h"

#include <stdio.h>
#include <stdlib.h>

void fatal_error(int const condition, char const * const msg, int const signal){
  if( condition ){
    printf("%s", msg);
    exit(signal);
  }
}
