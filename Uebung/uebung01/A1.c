/***********************************************************************
 *
 * kompilieren  mit
 *
 *   gcc -Wall -pedantic A1.c -o A1 -lm
 *
 *   Die Optionen fuer gcc koennen mit "man gcc" nachgelesen werden.
 *
 * run with
 *   ./A1
 *
 ***********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/***********************************************************************
 * Einlesen von dataset
 *
 * Array dout speichert die Daten, die aus Datei "filename" gelesen werden
 * we pass pointer to ( pointer to double ), so we can set the content
 *
 * dout und nout werden als Zeiger uebergeben ( Adresse von Zeiger dout und nout ),
 * um die Manipulationen des Inhalts auch ausserhalb der Funktion read_dataset
 * zu erhalten ( ansonsten Aenderungen der Funktionsargumente nur innerhalb der Funktion wirksam )
 *
 ***********************************************************************/
int read_dataset ( double ** const dout, unsigned int * const nout , char * const filename, int const binary ) {

  FILE * fp = fopen ( filename, "r" );
  if ( fp == NULL ) return ( 1 );

  unsigned int n = 0;
  double * d = NULL;

  /* Einlesen der Meta-Daten ( Daten zur Beschreibung des Datensatzes ) */
  if ( binary ) {
    if ( fread ( &n, sizeof ( unsigned int ), 1, fp ) !=  1 ) return( 2 );
  } else {
    if ( fscanf ( fp, "%u\n", &n ) != 1 ) return ( 2 );
  }
  /* Falls Null gesesen wurde, Abbruch */
  if ( n == 0 ) {
    fprintf ( stderr, "[read_dataset] no data\n" );
    return ( 4 );
  } else {
    fprintf ( stdout, "# [read_dataset] length of dataset = %d\n", n );
  }

  /* Array der Laenger n doubles anlegen */
  d = (double*) malloc ( n * sizeof(double) );
  if ( d == NULL ) {
    fprintf ( stderr, "[read_dataset] Error from malloc %s %d\n", __FILE__, __LINE__ );
  }
  /* Daten einlesen */ 
  if ( binary ) {
    /* Lesen im Binaerformat mit fread */
    if ( fread ( d, sizeof ( double ), n , fp ) !=  n ) return( 3 );
  } else {
    /* Lesen im ASCII-Format mit fscanf, zeilenweise */
    for ( unsigned int i = 0; i < n; i++ ) {
      if ( fscanf ( fp, "%lf\n", d+i ) != 1 )  return( 3 );
    }
  }  /* end of if binary */

  /* Schliesse Dateizeiger */
  fclose ( fp );

  /* Zuweisung der Rueckgabewerte */
  *nout = n;
  *dout = d;

  /* Rueckker aus der Funktion mit Rueckgabewert 0 = kein Fehler */
  return ( 0 );

}  /* end of read_dataset */

/***********************************************************************
 * Zeige den Inhalt des Feldes data in stdout
 ***********************************************************************/
int show_dataset ( double * const data, unsigned int const n ) {

  if ( data == NULL || n == 0 ) return ( 1 );

  fprintf( stdout, "\n\n# [show_dataset]\n" );
  for ( unsigned int i = 0; i < n; i++ ) {
    fprintf( stdout, "%6d %25.16e\n", i, data[i] );
  }
  return ( 0 );
}  /* end of show_dataset */

/***********************************************************************
 * kumulative Quadratsumme
 *
 * benutzt 2 Felder
 *   Eingangswerte in din
 *   Ausgangswerte in dout
 *
 * kann auch "in-place" genutzt werden, also dout = din und dout wird mit
 * den neu berechneten Werten ueberschrieben
 *
 ***********************************************************************/
int cum_sqr_sum ( double * const dout, double * const din, unsigned int n ) {
  if ( dout == NULL || din == NULL || n == 0 ) return ( 1 );

  dout[0] = din[0] * din[0];
  for ( unsigned int i = 1; i < n; i++ ) {
    dout[i] = dout[i-1] + din[i] * din[i];
  }

  return ( 0 );
} /* end of cum_sqr_sum */


/***********************************************************************
 * MAIN PROGRAM
 ***********************************************************************/

int main(int argc, char **argv) {

  unsigned int ndata = 0;
  int exitstatus;
  double * data = NULL;
  FILE * fout = NULL;

  /* dataset1 einlesen */
  exitstatus =read_dataset ( &data, &ndata, "dataset1", 0 );
  if ( exitstatus != 0 ) {
    fprintf( stderr, "[main] Error from read_dataset, status was %d %s %d\n", exitstatus, __FILE__, __LINE__);
    exit( 3 );
  }

  /* data in stdout anzeigen */
  show_dataset( data, ndata );

  /* kumulative Quadratsumme in data2 */
  double * data2 = (double*)malloc ( ndata * sizeof(double) );
  /* check, that array has been allocated */
  if (data2 == NULL ) exit( 5 );

  /* Datei oeffnen zum Schreiben */
  fout = fopen ( "dataset1.cum_sqr_sum", "w" );
  if ( fout == NULL) exit ( 6 );

  /* Funktionsaufruf */
  cum_sqr_sum ( data2, data, ndata );
  
  /* data2 nach fout schreiben */
  for( unsigned int i = 0; i < ndata; i++ ) {
    fprintf ( fout, "%e\n", data2[i] );
  }
  /* Datei schliessen */
  fclose ( fout );

  /* dynamisch allokierten heap Speicher wieder freigeben; (pruefe hier noch mal, dass Zeiger nicht NULL sind ) */
  if ( data  != NULL ) free ( data  );
  if ( data2 != NULL ) free ( data2 );

  /* das gleiche fuer dataset2 */
  exitstatus =read_dataset ( &data, &ndata, "dataset2", 1 );
  if ( exitstatus != 0 ) {
    fprintf( stderr, "[main] Error from read_dataset, status was %d %s %d\n", exitstatus, __FILE__, __LINE__);
    exit( 4 );
  }

  /* show the binary data set in stdandard output */
  show_dataset( data, ndata );

  /* apply sin function, save table */
  fout = fopen ( "dataset2.sin", "w" );
  if ( fout == NULL) exit ( 7 );

  /* Schleife ueber Feldelemente, wende Funktion sin() aus der Mathe-Bibliothek math an */
  for( unsigned int i = 0; i < ndata; i++ ) {
    fprintf ( fout, "%25.16e %25.16e\n", data[i], sin( data[i] ) );
  }
  fclose ( fout );

  /* heap-Speicher freigeben */
  if ( data != NULL ) free ( data );

  return ( 0 );
}
