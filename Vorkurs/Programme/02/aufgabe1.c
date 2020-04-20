#include <stdio.h>

int main(){
    int  n =  10;
    int  i;
    int summe =  0;
    
    i=0;

    while (i<=n){
        summe += i;
        i++;
    }
    printf("Das Ergebnis ist %d\n",  summe);
    return 0;
}