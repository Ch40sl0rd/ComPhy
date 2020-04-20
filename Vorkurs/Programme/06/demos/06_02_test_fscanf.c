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
  char string[100];

  int n_gelesen = 0;  

  char const * dateiname = "format.txt";
  FILE* fp = fopen(dateiname,"r");
  if( fp == NULL ){
    printf("\'%s\' konnte nicht geoeffnet werden!\n", dateiname);
    exit(21);
  }

  double x = 0.0;
 
  printf("Dateicursor bei: %ld\n", ftell(fp) );
  n_gelesen = fscanf(fp, "%s %lf\n", string, &x);
  printf("%s\n", string);
  printf("n_gelesen = %d\n", n_gelesen);
  printf("x = %lf\n",x);
  printf("Dateicursor bei: %ld\n", ftell(fp) );
  ende_erreicht_abfrage(fp); 

  x = 4.2;
  // wird fehlschlagen, falscher Format-string
  n_gelesen = fscanf(fp, "Test %lf\n", &x);
  printf("n_gelesen = %d\n", n_gelesen);
  printf("x = %f\n",x);
  printf("Dateicursor bei: %ld\n", ftell(fp) );
  ende_erreicht_abfrage(fp); 
  
  // versuchen wir den vermeintlich richtigen Format-string 
  n_gelesen = fscanf(fp, "Testb %lf\n", &x);
  printf("\nn_gelesen = %d\n", n_gelesen);
  printf("x = %f\n",x);
  printf("Dateicursor bei: %ld\n", ftell(fp) );
  ende_erreicht_abfrage(fp); 
  
  // ah, wir hatten 'Test' ja schon gelesen... 
  n_gelesen = fscanf(fp, "b %lf\n", &x);
  printf("n_gelesen = %d\n", n_gelesen);
  printf("x = %f\n",x);
  printf("Dateicursor bei: %ld\n", ftell(fp) );
  ende_erreicht_abfrage(fp); 

  int error = fclose(fp);
  if( error ){
    printf("Fehler beim Schliessen der Datei!\n");
    exit(1);
  }
  
  return 0;
}
