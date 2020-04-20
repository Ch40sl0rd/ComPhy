#include <stdio.h>
#include <unistd.h>
int main(void){
  char const * infilename = "sw.txt";
  char const * outfilename = "sw_copy.txt";

  FILE * infile = fopen(infilename, "r");
  FILE * outfile = fopen(outfilename, "w");

  if( (void*)infile == NULL ){
    printf("\'%s\' konnte nicht geoeffnet werden...\n", infilename);
    return 1;
  }
  if( (void*)outfile == NULL ){
    printf("\'%s\' konnte nicht geoeffnet werden...\n", outfilename);
    return 1;
  }

  int in = 0;
  while( (in = fgetc(infile)) != -1 ){
    putchar(in);
    if( (char)in == '\n' ){
      putchar('\n');
      sleep(1);
      putchar('\n');
    } else {
      fputc(in, outfile);
    }
  }

  // Noch explizit auf Fehler ueberpruefen!
  if( ferror(infile) ){
    printf("Es sind Fehler beim Lesen von \'%s\' aufgetreten!\n", infilename);
  }
  if( ferror(outfile) ){
    printf("Es sind Fehler beim Schreiben von \'%s\' aufgetreten!\n", outfilename);
  }

  int in_error = fclose(infile);
  int out_error = fclose(outfile);
  if( in_error ){
    printf("Fehler beim Schliessen der Eingabedatei!\n");
  }
  if( out_error ){
    printf("Fehler beim Schliessen der Ausgabedatei!\n");
  }
  if( in_error || out_error ){
    return(2);
  }

  return 0;
}
