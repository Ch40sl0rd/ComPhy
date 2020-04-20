#include <stdlib.h>
#include <stdio.h>
void ende_erreicht_abfrage(FILE * stream){
  printf("Ende der Datei erreicht? ");
  if( feof(stream) ){
    printf("Ja!\n\n");
  } else {
    printf("Nein!\n\n");
  }
}  

int main(void){
  char const * dateiname = "format.txt";
  FILE* fp = fopen(dateiname,"rb");
  if( fp == NULL ){
    printf("\'%s\' konnte nicht geoeffnet werden!\n", dateiname);
    exit(21);
  }

  double x = 0.0;
  int n_gelesen = 0;
 
  printf("Dateicursor bei: %ld\n", ftell(fp) );
  n_gelesen = fscanf(fp, "Test %lf\n", &x);
  printf("n_gelesen = %d\n", n_gelesen);
  printf("x = %lf\n",x);
  printf("Dateicursor bei: %ld\n", ftell(fp) );
  ende_erreicht_abfrage(fp); 

  // Dateicursor an den Anfang der Datei bewegen
  fseek(fp, 0, SEEK_SET);
  printf("Dateicursor bei: %ld\n", ftell(fp) );
  x = 4.2;

  // wir sind jetzt wieder am Anfang der Datei, also nochmal 'x' auslesen!  
  printf("Dateicursor bei: %ld\n", ftell(fp) );
  n_gelesen = fscanf(fp, "Test %lf\n", &x);
  printf("n_gelesen = %d\n", n_gelesen);
  printf("x = %lf\n",x);
  printf("Dateicursor bei: %ld\n", ftell(fp) );
  ende_erreicht_abfrage(fp); 
  
  // Jetzt versuchen wir, wie in Beispiel 06_02, einen falschen Format-string
  // zu nutzen
  printf("Dateicursor bei: %ld\n", ftell(fp) );
  n_gelesen = fscanf(fp, "Test %lf\n", &x);
  printf("n_gelesen = %d\n", n_gelesen);
  printf("x = %lf\n",x);
  printf("Dateicursor bei: %ld\n", ftell(fp) );
  ende_erreicht_abfrage(fp); 

  // Wir hängen jetzt an Byte 16 fest, könnten also 4 Bytes zurückhüpfen
  // und den Einleseversuch wiederholen
  // dazu bewegen wir uns 4 Bytes zurück von der momentanen Position 'SEEK_CUR'
  fseek(fp, -4, SEEK_CUR);
  printf("Dateicursor bei: %ld\n", ftell(fp) );
  
  printf("Dateicursor bei: %ld\n", ftell(fp) );
  n_gelesen = fscanf(fp, "Testb %lf\n", &x);
  printf("n_gelesen = %d\n", n_gelesen);
  printf("x = %lf\n",x);
  printf("Dateicursor bei: %ld\n", ftell(fp) );
  ende_erreicht_abfrage(fp); 
  

  // Jetzt begeben wir uns noch einen Byte vom Anfang der Datei 
  // via 'SEEK_SET'
  fseek(fp, 1, SEEK_SET);
  x = 4.2;
  printf("Dateicursor bei: %ld\n", ftell(fp) );
  // wir haben 'T' aus 'Test' übersprungen -> format string ist also
  // 'est %lf'
  n_gelesen = fscanf(fp, "est %lf\n", &x);
  printf("n_gelesen = %d\n", n_gelesen);
  printf("x = %lf\n",x);
  printf("Dateicursor bei: %ld\n", ftell(fp) );

  int error = fclose(fp);
  if( error ){
    printf("Fehler beim Schliessen der Datei!\n");
    exit(1);
  }
  
  return 0;
}
