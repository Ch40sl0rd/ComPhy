// die Praeprozessorkonstante '_POSIX_C_SOURCE' bestimmt, welcher "POSIX-Standard"
// von uns genutzt wird
#define _POSIX_C_SOURCE 200809L
// Die seit 200809 eingefuehrten Funktionen koennen jetzt auch genutzt werden
#include <stdio.h>

#include <stdlib.h>
#include <string.h>
int main(void){
  char * lineptr = NULL;
  size_t n_bufsize = 0;
  ssize_t n_read = 0;

  unsigned int const n_dateien = 3;
  unsigned int const max_laenge_dateiname = 100;

  char dateinamen[n_dateien][max_laenge_dateiname];

  snprintf(dateinamen[0], max_laenge_dateiname, "format.txt");
  snprintf(dateinamen[1], max_laenge_dateiname, "sw.txt");
  snprintf(dateinamen[2], max_laenge_dateiname, "pi_10000.txt"); 
  
  for( unsigned int i_datei = 0; i_datei < n_dateien; i_datei++ ){
    FILE *stream = fopen(dateinamen[i_datei], "r");
    if( stream == NULL ){
      printf("\'%s\' konnte nicht geoeffnet werden!\n", dateinamen[i_datei]);
      return 1;
    }

    unsigned int linecounter = 0;
    while( ( n_read = getline(&lineptr, &n_bufsize, stream) ) != -1 ){
      linecounter++;
      printf("Zeile %u: Es wurden %lu Zeichen gelesen\n", linecounter, n_read);
      printf("Zeile %u: Der Lesepuffer hat Groesse %lu\n", linecounter, n_bufsize);
      printf("Zeile %u: %s\n", linecounter, lineptr);
    }
    int error = fclose(stream);
    if( error ){
      printf("Fehler beim Schliesen von %s!\n", dateinamen[i_datei]);
    }

    printf("Taste druecken fuer naechste Datei!\n");
    getchar();
  }
  // zwischendurch kuemmert sich getline um den Speicher, am Ende muessen
  // wir aber den letzten Puffer freigeben
  free(lineptr);

  return 0;
}
