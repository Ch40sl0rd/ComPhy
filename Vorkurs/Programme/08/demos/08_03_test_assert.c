#include <assert.h>
#include <stdio.h>
int func(FILE* fp){
  assert( fp != NULL );
  /* ... */
  return 0;
}
int main(void){
  FILE* fp = fopen("gibtbesnicht.txt","r");
  func(fp);
  return 0;
}
