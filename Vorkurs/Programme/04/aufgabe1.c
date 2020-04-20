#include<stdio.h>
void tausch( int* a, int* b){
    //Hier werden zwei Pointer benötigt und deren Inhalt vertauscht
    int zwischen = (*a);
    //Wir können pointer derefenzieren mit(*a), damit können wir auf die Speicherstelle, auf die der Pointer zeigt, zugreifen.
    (*a) = (*b);
    (*b) = zwischen;
}

int main(){
    int a = 5;
    int b = 10;
    printf("Die Werte von a und b sind: %d und %d.\n", a, b);
    tausch(&a, &b);
    printf("Die Werte von a und b sind: %d und %d.\n", a, b);
    return 0;
}