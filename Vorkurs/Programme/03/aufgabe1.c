#include <stdio.h>

long long fakultaet(int n){
    long  long ergebnis = 1;
    while(n>1){
        ergebnis*=n;
        n--;
    }
    return ergebnis;
}

int main(){
    long long add2fak;
    add2fak = fakultaet(5) + fakultaet(17);
    printf("5! + 17! = %lld\n", add2fak);
    return 0;
}