#ifndef POS_ST_H
#define POS_ST_H

#include <math.h>

typedef struct pos_t {
  double x;
  double y;
  double z;
} pos_t;

// Distanz zum Urprung.
double pos_urspr_dist(pos_t const * const pos);

// Gibt zurueck, ob 'links' oder 'rechts' naeher am Urpsrung liegt. 
// Entspricht der Spezifikation fuer qsort aus
// der C Standardbibliothek.
int pos_compar_urspr_dist( void const * links, void const * rechts );

// Gibt zurueck, ob 'links' oder 'rechts naeher an der XY-Ebene liegt
// Entspricht der Spezifikation fuer qsort aus
// der C Standardbibliothek.
int pos_compar_xy_ebene_dist( void const * links, void const * rechts );

// Initialisiert ein Array aus N pos_t zufaellig im 3D Kubus [-1,-1,-1] bis [1,1,1]
void pos_init_rand( pos_t * posarray, unsigned int const N);

// Gibt Positionskoordinaten sowie Distanz zum Ursprung auf dem Bildschirm aus
void pos_printf( pos_t const * const pos );

#endif // POS_ST_H
