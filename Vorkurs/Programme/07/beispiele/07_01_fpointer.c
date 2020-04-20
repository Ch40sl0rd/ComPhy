#include <stdio.h>

double func(double omega){
  return 2*omega;
}

double void_func(void){
  return (double)3234.32;
}

int main(void){
  // wir werden versuchen, zu verstehen, was schiefgehen kann, wenn man
  // die falschen Funktionenzeiger nutzt. Die Ergebnisse dieses Tests
  // sind aller Wahrscheinlichkeit nach systemabhängig
  double (*void_fp)(void);
  double (*fp)(double);

  // Was passiert, wenn wir die falschen Funktionenpointer nutzen?
  // Compiler wird warnen, aber keinen Fehler produzieren..
  fp = void_func;
  void_fp = func;

  double w, x, y, z1, z2;

  // void_func direkt oder via fp aufrufen
  w = void_func();
  x = (*fp)(32); // Argument wird ignoriert werden

  printf("w = %f\n", w);
  printf("x = %f\n", x);
 
  // was passiert, wenn wir die Funktion 'func' über einen 
  // Funktionenpointer ohne Argumente aufrufen? 
  z1 = (*void_fp)();
  printf("erster Aufruf von func ueber void_fp: z1 = %f\n", z1);

  // rufen wir mal func direkt auf, aber ohne, dass wir auf den
  // Rückgabewert zugreifen
  func(24);

  // und jetzt rufen wir 'func' noch ein zweites Mal über 'void_fp' auf
  z1 = (*void_fp)();

  // z1 wird jetzt 96.0 sein! Der sogenannte 'Stackpointer' liegt zufällig
  // Rückgabewert des Aufrufs von 'func' mit dem Argument '24' in Zeile 36
  // Dieser Wert wird dan von 'func' als Wert für 'omega' genutzt, wenn 'func'
  // in Zweile 37 über 'void_fp' ohne Argumente aufgerufen wird
  printf("zweiter Aufruf von func ueber void_fp: z1 = %f\n", z1);

  // versuchen wir das nochmal, aber diesmal mit Zuweisung
  z1 = func(24);
  printf("zweiter direkter Aufruf von func: z1 = %f\n", z1);

  // und jetzt rufen wir 'func' noch ein drittes Mal über 'void_fp' auf
  // in meinem Fall ist 'z1' dann 0.0, da wir in Zeile 51 den Rückgabewert
  // von func ausgelesen haben
  z1 = (*void_fp)();
  printf("dritter Aufruf von func ueber void_fp: z1 = %f\n", z1);

  // Zum Glück gibt es hier zumindest eine Fehlermeldung, da wir hier
  // versuchen über einen void Funktionenpointer ein Argument
  // zu übergeben
  //z1 = (*void_fp)(55);

  fp = func;
  void_fp = void_func;

  y = (*fp)(32); // ah, jetzt geht es!
  z2 = (*void_fp)();

  printf("y =  %f, z1 =  %f, z2 = %f\n", y, z1, z2);

  return 0;
}
