#include <stdlib.h>
#include <stdio.h>
int main(void){
  const unsigned int laenge = 40;

  double *x = malloc( laenge*sizeof(double) );
  if( (void*)x == NULL ){
    printf("Speicherallokation fehlgeschlagen!\n");
    return 1;
  }

  char const * dateiname = "x.dat";
  char const * modus = "rb";

  FILE* eingabedatei = fopen(dateiname, modus);
  if( eingabedatei == NULL ){
    printf("\'%s\' konnte nicht im Modus \'%s\' geoeffnet werden!\n", dateiname, modus);
    return 1;
  }
  
  // hier muss man natuerlich wieder aufpassen: x muss gross genug sein,
  // sonst wird einfach unerlaubt in irgendwelchen Speicher geschrieben...
  unsigned int rval = fread( (void*)x, sizeof(double), laenge, eingabedatei );
  if( rval != laenge ){
    printf("Fehler beim Lesen\n");
  }
  if( ferror(eingabedatei) ){
    printf("ferror != 0\n");
  }
  printf("feof = %d\n", feof(eingabedatei) );
  int error = fclose(eingabedatei);
  if( error ){
    printf("Fehler beim Schliessen der Datei!\n");
    exit(1);
  }
  
  for(unsigned int i = 0; i < laenge; ++i ){
    printf("%u %f\n", i, x[i]);
  } 
  
  free(x);
  return 0;
}
