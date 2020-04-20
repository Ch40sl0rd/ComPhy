#include "teilchen_2d.h"

#include <stdlib.h>
#include <stdio.h>

int main(void){
  unsigned int const laenge = 4;

  teilchen_2d_t teilchen[laenge];

  char const * dateiname = "teilchen.dat";
  char const * modus = "rb";

  FILE* eingabedatei = fopen(dateiname, modus);
  if( eingabedatei == NULL ){
    printf("\'%s\' konnte nicht im Modus \'%s\' geoeffnet werden!\n", dateiname, modus);
    return 1;
  }

  unsigned int rval = fread( (void*)teilchen, sizeof(teilchen_2d_t), laenge, eingabedatei );
  if( rval != laenge ){
    printf("Fehler beim Lesen!\n");
  }
  if( ferror(eingabedatei) ){
    printf("ferror != 0\n");
  }
  printf("feof = %d\n", feof(eingabedatei) );
  int error = fclose(eingabedatei);
  if( error ){
    printf("Fehler beim Schliessen der Datei!\n");
    return 2;
  }
  
  teilchen_2d_print( &teilchen[2] );
  
  return 0;
}
