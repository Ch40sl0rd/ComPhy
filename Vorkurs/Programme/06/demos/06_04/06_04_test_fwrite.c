#include "teilchen_2d.h"

#include <stdlib.h>
#include <stdio.h>

int main(void){
  unsigned int const laenge = 4;

  teilchen_2d_t teilchen[laenge];
  teilchen_2d_init(teilchen, laenge);

  teilchen_2d_print(&teilchen[2]);

  char const * dateiname = "teilchen.dat";
  char const * modus = "wb";

  FILE * ausgabedatei = fopen(dateiname, modus);
  if( ausgabedatei == NULL ){
    printf("\'%s\' konnte nicht im Modus \'%s\' geoeffnet werden!\n", dateiname, modus);
    return 1;
  }

  unsigned int rval = fwrite( (void*)teilchen, sizeof(teilchen_2d_t), laenge, ausgabedatei );
  if( rval != laenge ){
    printf("Fehler beim Schreiben\n");
  }
  if( ferror(ausgabedatei) ){
    printf("ferror != 0\n");
  }
  printf("feof = %d\n", feof(ausgabedatei) );
  int error = fclose(ausgabedatei);
  if( error ){
    printf("Fehler beim Schliessen der Datei!\n");
    return 2;
  }

  return 0;
}
