#ifndef DLL_H
#define DLL_H
#include <stdlib.h>
#include <stdio.h>

typedef struct dll_node_t{
    struct dll_node_t *prev;
    struct dll_node_t *next;
    double  data;
}dll_node_t;

typedef struct dll_t{
    struct dll_node_t *first;
    struct dll_node_t *last;
}dll_t;

/* 
    Diese Funktion erzeugt eine neue, leere Liste.

    @list (out): list-pointer.

    funktioniert.
*/
dll_t*  dll_create(){
    dll_t* list = (dll_t*)malloc(sizeof(dll_t));
    if(!list){
        exit(1);
    }
    list->first = NULL;
    list->last = NULL;
    return list;
}

/*
    Diese Funktion löscht eine gesamte Liste und ihren Inhalt.

    @list (in): Liste, die gelöscht werden soll.

    funktioniert.
*/
void dll_free(dll_t *list){
    if(list==NULL){
        printf("Es existiret keine Liste.\n");
        return;
    }
    if(list->first == NULL){
        free(list);
        printf("Die Liste ist bereits leer und kann nicht gelöscht werden.\n");
        return;
    }
    dll_node_t *temp = list->first;
    dll_node_t *temp2;
    while(temp != NULL){
        temp2 = temp->next;
        free(temp);
        temp = temp2;
    }
    free(list);
}

/*
    Diese Funktion fügt einen neuen Wert nach dem gegeben Knoten ein.
    Dabei bedeutet der Knoten NULL den Anfang der Liste.

    @list (in): Liste zu der etwas hinzugefügt werden soll.
    @cursor (in): Knoten, nach dem der neue Wert eingefügt werden soll.
    @data (in): Wert, der eingefügt werden soll.

    @list_node (out): Der neue Knoten, der eingefügt wurde.

    funktioniert.
*/
dll_node_t* dll_insert(dll_t *list, dll_node_t *cursor, double data){
    if(list == NULL){
        printf("Es wurde keine Liste gefunden.\n");
        exit(2);
    }
    dll_node_t *neu = (dll_node_t*)malloc(sizeof(dll_node_t)*1);
    if(!neu){
        printf("Fehler bei Speicherallokierung!\n");
        exit(1);
    }
    neu->data = data;
    neu->next = NULL;
    neu->prev = NULL;

    if(cursor == NULL){
        if(list->last == NULL){
            list->first = neu;
            list->last = neu;
        }
        else{
            neu->next = list->first;
            list->first->prev = neu;
            list->first =  neu;
        }
        return neu;
    }
    
    if(cursor->next == NULL){
        cursor->next = neu;
        neu->prev = cursor;
        list->last = neu;
        return neu;
    }

    dll_node_t *temp = cursor->next;
    cursor->next = neu;
    neu->prev = cursor;
    neu->next = temp;
    return neu;
}

/*
    Diese Funktion löscht den angegebenen Knoten aus der gegebenen Liste.

    @list (in): Liste, deren Knoten gelöscht werden soll.
    @delnote (in): Knoten, der aus der Liste entfernt werden soll.

    funktioniert soweit.
*/
void dll_delete(dll_t *list, dll_node_t *delnote){
    if(list==NULL){
        printf("Es konnte keine Liste gefunden werden!\n");
        exit(2);
    }
    if(delnote == list->first){
        if(delnote == list->last){
            free(delnote);
            list->first = NULL;
            list->last = NULL;
        }
        delnote->next->prev = NULL;
        list->first = delnote->next;
        free(delnote);
        return;
    }
    if(delnote ==  list->last){
        delnote->prev->next =  NULL;
        list->last =  delnote->prev;
        free(delnote);
        return;
    }
    dll_node_t *temp  = delnote;
    delnote->next->prev = temp->prev;
    delnote->prev->next = temp->next;
    free(delnote);
    
}

/*
    Diese Funktion fügt einen Wert an den Anfang der Liste ein.

    @list (in): liste, an deren Anfang geschrieben werden soll.
    @data (in): Zahl, die an den Anfagn geschrieben werden soll.

    @list_node (out): Neuer Datenknoten, der nun Am Ende der Liste steht

    soweit getestet.
*/ 
dll_node_t* dll_unshift(dll_t *list, double data){
    if(list==NULL){
        printf("Keine Liste gefunden!\n");
        exit(404);
    }
    dll_node_t *neu = (dll_node_t*)malloc(sizeof(dll_node_t));
    if(!neu){
        printf("Fehler bei Speicherallokierung!\n");
        exit(1);
    }
    neu->data = data;
    neu->prev = NULL;
    neu->next = NULL;

    if(list->first == NULL){
        printf("Das ist eine leere Liste.\n");
        list->first = neu;
        list->last = neu;
        return neu;
    }
    neu->next = list->first;
    list->first->prev = neu;
    list->first = neu;
    return neu;
}
/*
    Diese Funktion fügt einen neuen Wert ans Ende der Liste ein.

    @list (in): Liste an deren Ende geschrieben werden soll.
    @data (in): Daten, die ans Ende geschrieben werden sollen.

    @list_node (out): Knoten der nun Am Ende der Liste steht.

    soweit getestet.
*/
dll_node_t* dll_push(dll_t *list, double data){
    if(list == NULL){
        printf("Keine Liste gefunden!\n");
        exit(404);
    }
    dll_node_t *neu = (dll_node_t*)malloc(sizeof(dll_node_t));
    if(!neu){
        printf("Fehler bei Speicherallokierung!\n");
        exit(1);
    }
    //Initalisierung des neuen Knoten
    neu->data = data;
    neu->prev = NULL;
    neu->next = NULL;
    //Auf leere Liste überprüfen
    if(list->first  == NULL){
        //printf("Es handelt sich um eine leere Liste!\n");
        list->first = neu;
        list->last = neu;
        return neu;
    }
    //Sonst falls es keine leere Liste ist.
    neu->prev = list->last;
    list->last->next = neu;
    list->last = neu;
    return neu;
}

/*
    Diese Funktion entfernt den Wert an erster Stelle der Liste und
    gibt diesen zurück.

    @list (in): liste, deren erstes Element gelöscht werden soll.

    @data  (out): Wert des ersten Elements;

    funktioniert soweit
*/
double dll_shift(dll_t *list){
    if(list->first == NULL){
        printf("Die Liste ist leer.\n");
        return -1.;
    }
    dll_node_t *cur = list->first;
    double ret = cur->data;
    //Prüfe ob die Liste nur ein Element enthält.
    if(cur->next == NULL){
        printf("Jetzt ist die Liste leer.\n");
        list->first = NULL;
        list->last = NULL;
        free(cur);
        return ret;
    }
    //Sonst entferne nur das erste Elment.
    list->first->next->prev = NULL;
    list->first = cur->next;
    free(cur);
    return ret;
}

/*
    Daten am Ende der Liste entfernen und zurückgeben.
    @(dll_t *list)in : double-linked-list pointer

    @data (out): Wert am Ende der Liste.

    fuunktioniert soweit.
*/
double dll_pop(dll_t *list){
    if(list == NULL){
        printf("Es konnte keine Liste gefunden werden.\n");
        exit(404);
    }
    dll_node_t *cur = list->last;
    
    if(!cur){
        printf("Die Liste ist leer.\n");
        return -1.;
    }

    double ret = cur->data;
    if(cur->prev == NULL){//Die Liste enthält dann nur ein Element
        list->first = NULL;
        list->last = NULL;
        free(cur);
        return ret;
    }
    cur->prev->next = NULL;
    list->last = cur->prev;
    free(cur);
    return ret;
}

/*
    Merge von zwei Listen, wobei die zweite Liste an die erste angefügt wird
    Die beiden Listpointer werden dabei zerstört und nur eine neue gemeinsame
    Liste zurückgegeben.

    @list1 (in): erste Liste, wird zu anfang gestellt.
    @list2 (in): zweite Liste, wird an erste Liste angefügt.

    @out: list-pointer auf gemeinsame Liste.
*/
dll_t* dll_merge(dll_t *list1, dll_t *list2){
    if( (list1 == NULL) || (list1->first == NULL) ){
        printf("Die erste Liste ist leer.\n");
        free(list1);
        return list2;
    }
    if((list2 == NULL) || (list2->first==NULL)){
        printf("Die zweite Liste ist leer.\n");
        free(list2);
        return list1;
    }
    dll_t* new_list = (dll_t*)malloc(sizeof(dll_t));
    if(!new_list){
        printf("Fehler bei Speicherallokierung!\n");
        exit(1);
    }
    list1->last->next = list2->first;
    list2->first->prev = list1->last;
    new_list->first = list1->first;
    new_list->last = list2->last;
    return new_list;
}

/*
    Ausgabe eine Liste,entweder vor- oder rückwärts.

    @list (in): liste die ausgegeben werden  soll.
    @option (in): vor(1)- oder rück(-1)-wärts, Richtung der Ausgabe.

    funktioniert momentan nur vorwärts.
*/
void dll_print(dll_t *list, int option){
    if(list == NULL){
        printf("Keine Liste gefunden!\n");
        exit(404);
    }
    if( (option*option)  > 1){
        printf("Keine gültige Option.\n");
        return;
    }
    if(option == 1){
        dll_node_t *current;
        if(list->first == NULL){
            printf("Die Liste ist leer.\n");
            return;
        }
        current = list->first;
        while(current->next != NULL){
            printf("%f; ", current->data);
            current = current->next;
        }
        printf("%f\n", current->data);
    }
}
#endif