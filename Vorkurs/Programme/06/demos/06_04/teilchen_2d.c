#include "teilchen_2d.h"

#include <math.h>
#include <stdio.h>

void 
teilchen_2d_init(teilchen_2d_t * teilchen,
                 const unsigned int N)
{
  for( unsigned int i = 0; i < N; ++i ){
    teilchen[i].x      = cos((double)i / 15);
    teilchen[i].y      = sin((double)i / 4);
    teilchen[i].v_x    = cos((double)i / 22) + sin((double)i / 4);
    teilchen[i].v_y    = sin((double)i / 11);
    teilchen[i].ladung = -1.0;
    teilchen[i].m      = 12.3;
  }
}

void
teilchen_2d_print(const teilchen_2d_t * teilchen){
  printf("  x = %.3e,   y = %.3e\n"
         "v_x = %.3e, v_y = %.3e\n"
         "     m = %.3e\n"
         "ladung = %.3e\n",
         teilchen->x,
         teilchen->y,
         teilchen->v_x,
         teilchen->v_y,
         teilchen->m,
         teilchen->ladung);
}
