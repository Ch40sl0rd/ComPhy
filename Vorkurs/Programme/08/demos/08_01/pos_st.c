#include "pos_st.h"

#include <stdlib.h>
#include <stdio.h>

double pos_urspr_dist( pos_t const * const pos ){
  return sqrt( pos->x*pos->x + pos->y*pos->y + pos->z*pos->z );
}

int pos_compar_urspr_dist( void const *links, void const *rechts ){
  double d_links = pos_urspr_dist( (pos_t const*)links );
  double d_rechts = pos_urspr_dist( (pos_t const*)rechts );
  if( d_links > d_rechts ){
    return 1;
  } else if ( d_links < d_rechts ){
    return -1;
  } else {
    return 0;
  }
}

int pos_compar_xy_ebene_dist( void const *links, void const *rechts ){
  double d_links = fabs( ((pos_t const *)links)->z );
  double d_rechts = fabs( ((pos_t const *)rechts)->z );
  if( d_links > d_rechts ){
    return 1;
  } else if ( d_links < d_rechts ){
    return -1;
  } else {
    return 0;
  }
}

void pos_init_rand( pos_t * poslist, unsigned int const N){
  // wir nutzen die Zufallszahlen hier ohne den Seed festzulegen
  // in der Praxis sollte man den Zufallszahlengenerator aber 
  // natuerlich initialisieren!
  // rand() ist ein sehr schlecher Zufallszahlengenerator, in der Praxis
  // lieber einen Generator aus GSL nutzen oder alternativ 
  // Ranlux oder Mersenne Twister direkt nehmen!
  for( unsigned int i = 0; i < N; ++i ){
    // zufaellige Koordinaten im 3D Kubus [-1,-1,-1] bis [1,1,1]
    poslist[i].x = -1.0 + 2*((double)rand())/RAND_MAX;
    poslist[i].y = -1.0 + 2*((double)rand())/RAND_MAX;
    poslist[i].z = -1.0 + 2*((double)rand())/RAND_MAX;
  }
}

void pos_printf( pos_t const * pos ){
  printf("x = %f, y = %f, z = %f\n", pos->x, pos->y, pos->z);
  printf("|r| = %f\n", pos_urspr_dist( pos ) );
}

