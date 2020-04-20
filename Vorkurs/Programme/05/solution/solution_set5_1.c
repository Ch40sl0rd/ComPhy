// (C) Carsten Urbach, 2017

#include<stdlib.h>
#include<stdio.h>

// wir benutzen explizit int!
// mit long unsigned int funktionert es natuerlich nicht!
typedef int myint;

// Diese Funktion erzeugt die Zufallszahlen
// Eingabe ist die vorige Zufallszahl
// Rueckgabe ist die naechste.
unsigned long int prng(const myint d) {
  const myint m = 2147483647;
  const myint a = 16807;
  const myint q = 127773;
  const myint r = 2836;

  myint res = a*(d % q) - r*(d/q);
  if(res < 0) {
    res += m;
  }
  return res;
}

int main() {
  // Anzahl an zu ziehenden Zufallszahlen
  int N = 10000;
  // erste Zufallszahl
  myint d = 232556;

  for(unsigned int i = 0; i < N; i++) {
    d = prng(d);
    printf("%d %f\n", d, (double)d/(double)2147483647);
  }
  return 0;
}
