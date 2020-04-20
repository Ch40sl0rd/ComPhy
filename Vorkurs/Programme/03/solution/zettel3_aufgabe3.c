#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <sys/time.h>
double power(double x, int n){
  if (n<0){
    printf("Fehler argument n soll >0 sein\n");
    exit(1);
  }
  if (n==0)
    return 1.;
  else if (n%2==0)
    return power(x*x, n/2);
  else 
    return x*power(x*x,(n-1)/2);
}
//Wir müssen extra speicher verwenden 
//ret-->fur das Ergebnis
//xi -->für die Aktuelle basis
//Beispiel 3**5
//ret=1 xi=3 i=5
//Esten Schritt -->>Wir mussen das ergebnis mit 3 multiplizieren,und wir mussen 9**2 berechnen: 
//ret=3 xi=9 i=2
//Zweite Schritt -->>Wir mussen 9**2 berechnen, laut double and add
//ret=3 xi=81 i=1 
//Dritte Schritt
//ret=243 xi=81**2 i=0 -->>Wir sind fertig
double xpower(double x, int n){
  double ret=1;
  double xi=x;
  int i=0;
  for (i=n;i>0;){
    if (i%2==0){
      xi=xi*xi;
      i=i/2;
    }
    else{
      ret=xi*ret;
      xi=xi*xi;
      i=(i-1)/2;
    }
  }
  return ret;
}
double naiv_power(double x,int n){
  double ret=1.;
  int i;
  for (i=0;i<n;++i)
   ret*=x;
  return ret;
}
int main(int argc,char *argv[]){
  struct timeval anfang, ende;
  gettimeofday(&anfang, NULL);
  double xnaiv;
  double x=xpower(0.9999999999,2000000000);
  gettimeofday(&ende, NULL);
  double sec = (double)(ende.tv_sec-anfang.tv_sec);
  double usec = (double)(ende.tv_usec-anfang.tv_usec);
  printf("Es sind %lf Sekunden vergangen double add\n", sec+1.0e-6*usec);
//  printf("Es sind %lf Sekunden vergangen in recursive\n", (double)(ende.tv_sec-anfang.tv_sec));

  gettimeofday(&anfang, NULL);
  xnaiv=power(0.9999999999,2000000000);
  gettimeofday(&ende, NULL);
  sec = (double)(ende.tv_sec-anfang.tv_sec);
  usec = (double)(ende.tv_usec-anfang.tv_usec);
  printf("Es sind %lf Sekunden vergangen naive\n", sec+1.0e-6*usec);
//  printf("Es sind %lf Sekunden vergangen in naive\n", (double)(ende.tv_sec-anfang.tv_sec));


  printf("Double add %e\n",x);
  printf("Naiv %e\n",xnaiv);
}
