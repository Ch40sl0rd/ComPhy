#include <stdlib.h>
#include <stdio.h>
int main(void){
  const unsigned int laenge = 40;

  double *x = malloc( laenge*sizeof(double) );
  for( unsigned int i = 0; i < laenge; i++ ){
    x[i] = (double)i*i*i;
  }

  char const * dateiname = "x.dat";
  char const * modus = "wb";

  FILE* ausgabedatei = fopen(dateiname, modus);
  if( ausgabedatei == NULL ){
    printf("\'%s\' konnte nicht im Modus \'%s\' geoeffnet werden!\n", dateiname, modus);
    return 1;
  }

  unsigned int rval = fwrite( (void*)x, sizeof(double), laenge, ausgabedatei );
  if( rval != laenge ){
    printf("Fehler beim Schreiben\n");
  }
  if( ferror(ausgabedatei) ){
    printf("ferror != 0\n");
  }
  printf("feof = %d\n", feof(ausgabedatei) );
  
  free(x);
  
  int error = fclose(ausgabedatei);
  if( error ){
    printf("Fehler beim Schliessen der Datei!\n");
    return 2;
  }
  
  return 0;
}
