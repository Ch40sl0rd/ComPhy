// B. Kostrzewa, M. Ueding, F. Pittler

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int str_compar(const void *links, const void *rechts) {
  // Bei den zu sortierenden Daten handelt es sich um Strings, also Arrays
  // von Zeichen.
  // 
  // Wenn wir jetzt einen Zeiger auf einen String haben, dann hat dieser Zeiger
  // den Datentyp 'char**'
  //
  // Die zu sortierenden Daten sind Arrays von Strings
  //
  // char * array[3];
  //
  // z.B., und haben auch Datentyp 'char**'.
  //
  // qsort liefert der 'compar' Funktion Zeiger auf zwei Elemente (links und rechts),
  // sagen wir mal das erste und zweite. Intern passiert folgendes in qsort:
  //
  // links = (const void*)(&array[0])
  // rechts = (const void*)(&array[1])
  //
  // also zwei Zeiger auf Strings mit zusaetzlichem 'const'
  //
  // Jetzt hat der void-Zeiger nur ein Sternchen, da es sich hierbei ja nur um einen
  // "beliebigen" Zeiger handelt. Wir müssen selbst wissen, wie wir diesen Zeiger
  // interpretieren müssen. (in der Vorlesung waren dies Zeiger auf structs)
  //
  // Hier ist die Interpretation: es handelt sich um einen Zeiger auf ein Element
  // aus einem Array von Strings. Also: char**
  //
  // ein wenig kompliziert wird dies durch die ganzen 'const'

  // Typecasten: 'links_cast' hat den Datentyp: ein Zeiger (nicht const) auf einen
  // Zeiger (auch nicht const), welcher selbst auf const-Daten zeigt (den Anfang eines Strings)
  char * const *links_cast = links;
  char * const *rechts_cast = rechts;

  // strcmp hat als Argumente 'const char*', also Strings beziehungsweise char-Zeiger
  // auf das erste Element eines konstanten Strings -> const char*
  // wir derferenziegen also die gecasteten Zeiger von oben und erhalten so
  // die eigentlichen Strings (also Zeiger auf das erste Element in einem char Array)
  const char *links_string = *links_cast;
  const char *rechts_string = *rechts_cast;

  // jetzt koennen wir strcmp aufrufen
  return strcmp(links_string, rechts_string);
}

int main() {

  const unsigned int N = 8;
  const unsigned int max_len =100;
  
  // Speicherplatz fuer 8 Zeiger auf Strings
  char **array = malloc(sizeof(char*)*N);
  if( (void*)array == NULL ){
    printf("Speicherallokation fuer 'char **array' fehlgeschagen!\n");
    return 1;
  }

  // Speicherplatz fuer die jeweiligen Strings
  for( unsigned int i = 0; i < N; ++i ){
    array[i] = malloc(sizeof(char)*max_len);
    if( (void*)array[i] == NULL ){
      printf("Speicherallokation fuer 'char*' (array[%d]) fehlgeschlagen!\n", i);
    }
  }
  
  // wir befuellen unsere 8 strings mit den Textschnipseln aus der Aufgabe
  snprintf(array[4], max_len, "%s", "1");
  snprintf(array[6], max_len, "%s", "10");
  snprintf(array[3], max_len, "%s", "10");
  snprintf(array[2], max_len, "%s", "123");
  snprintf(array[7], max_len, "%s", "2344");
  snprintf(array[0], max_len, "%s", "Hallo");
  snprintf(array[1], max_len, "%s", "Thor");
  snprintf(array[5], max_len, "%s", "Tor");

  for (unsigned int i = 0; i < N; ++i) {
    printf("%i: %s\n", i, array[i]);
  }
  printf("\n");

  qsort(array, N, sizeof(char*), str_compar);

  for (unsigned int i = 0; i < N; ++i) {
    printf("%i: %s\n", i, array[i]);
  }
  printf("\n");
  
  // Speicher fuer die einzelnen Strings freigeben
  for( unsigned int i = 0; i < N; ++i ){
    free(array[i]);
  }
  // jetzt erst Speicher fuer das Array von Strings freigeben
  free(array);

  return 0;
}

