#include "pos_st.h"

#include <stdio.h>
#include <stdlib.h>

int main(void){
  const unsigned int N = 20;
  pos_t positionen[N];
  pos_init_rand( positionen, N );

  for(unsigned int i = 0; i < N; ++i ){
    pos_printf( &positionen[i] );
  }

  printf("Positionen werden jetzt nach Distanz zum Ursprung sortiert!\nTaste druecken, um weiterzumachen.\n");
  getchar();

  qsort( (void*)positionen, 
         (size_t)N,
         sizeof(pos_t), 
         pos_compar_urspr_dist );

  for(unsigned int i = 0; i < N; ++i ){
    pos_printf( &positionen[i] );
  }
  
  printf("Positionen werden jetzt nach Distanz zur XY-Ebene  sortiert!\nTaste druecken, um weiterzumachen.\n");
  getchar();
  qsort( (void*)positionen, 
         (size_t)N,
         sizeof(pos_t), 
         pos_compar_xy_ebene_dist );

  for(unsigned int i = 0; i < N; ++i ){
    pos_printf( &positionen[i] );
  }

  return 0;
}

