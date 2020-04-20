#include <string.h>
#include <stdlib.h>
#include <stdio.h>

int main(void){
  const int n_pfade = 3;
  const int laenge = 200;
  char pfade[n_pfade][laenge];
  
  // Datei im momentanen Verzeichnis (pwd)
  snprintf(pfade[0], laenge, "%s", "sw.txt");
  // relativer Pfad im UNIX Format
  snprintf(pfade[1], laenge, "../06/%s", "sw.txt");
  // relativer Pfad im Windows Format
  // '\' ist ein spezielles Zeichen, also muss man es 'escapen', mit '\' -> '\\'
  snprintf(pfade[2], laenge, "..\\06\\%s", "sw.txt");

  FILE *sw_f;
  for( unsigned int i = 0; i < n_pfade; ++i ){
    sw_f = fopen( pfade[i], "r");
    if( sw_f == NULL ){
      printf("Oeffnen von \'%s\' fehlgeschlagen!\n", pfade[i]);
    } else {
      printf("Oeffnen von \'%s\' gelungen! Wird geschlossen.\n", pfade[i]);
      int error = fclose(sw_f);
      if( error ){
        printf("Schliessen von \'%s\' fehlgeschlagen!\n", pfade[i]);
        printf("Programm wird beendet!\n");
        return 1;
      }
    }
  }
  
  return 0;
}
