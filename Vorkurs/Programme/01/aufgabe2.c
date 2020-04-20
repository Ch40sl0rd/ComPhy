/* Hello World Programm
* (c) someone */

#include <stdio.h>
#include <math.h>

int func(int n){
    if(n%2==0) {
        return n/2;
    }
    else if (n%2==1)
    {
        return (n+1) / 2;
    }
    else{
        return -1;
    }
    
}

int alg1(int c){
    int n = 2;
    int rootc = (int)sqrt( (double)c);
    while(1){
        if(n>rootc+1){
            return 1;
        }
        else if(c%n==0){
            return 0;
        }
        n++;
    }
    return 0;
}

int main () {
    printf("Hallo Welt\n");
    int m = 1003;
    printf("Und der Wert der Funktion bei m=%d ist %d\n", m, func(m));
    printf("Wir testen ob m=%d eine Primzahl ist: %d\n", m, alg1(m));
    return 0;
}