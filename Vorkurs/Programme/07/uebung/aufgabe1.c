#include <stdlib.h>
#include <stdio.h>
#include "dll.h"


int main(){
    dll_t *liste = dll_create();
    dll_push(liste, 1.);
    dll_unshift(liste, 2.);
    dll_unshift(liste, 6.);
    dll_unshift(liste, 8.);
    dll_unshift(liste, 13.);
    dll_push(liste, 3.);
    dll_insert(liste, liste->first->next, 42.5);
    dll_insert(liste, liste->last, 9000.5);
    dll_print(liste, 1);

    dll_t *liste2 = dll_create();
    dll_push(liste2, 123.);
    dll_unshift(liste2, 321.);
    dll_print(liste2, 1);

    dll_t *newlist  = dll_merge(liste, liste2);
    dll_print(newlist, 1);

    
    printf("Der Wert an erster Stelle beträgt:  %f\n",  dll_shift(liste));
    printf("Der Wert an letzter Stelle beträgt:  %f\n",  dll_pop(liste));
    dll_print(liste, 1);
    dll_free(liste);
}