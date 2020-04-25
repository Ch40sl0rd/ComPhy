#ifndef DLL_H
#define DLL_H
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
dll_t*  dll_create();

/*
    Diese Funktion löscht eine gesamte Liste und ihren Inhalt.

    @list (in): Liste, die gelöscht werden soll.

    funktioniert.
*/
void dll_free(dll_t *list);

/*
    Diese Funktion fügt einen neuen Wert nach dem gegeben Knoten ein.
    Dabei bedeutet der Knoten NULL den Anfang der Liste.

    @list (in): Liste zu der etwas hinzugefügt werden soll.
    @cursor (in): Knoten, nach dem der neue Wert eingefügt werden soll.
    @data (in): Wert, der eingefügt werden soll.

    @list_node (out): Der neue Knoten, der eingefügt wurde.

    funktioniert.
*/
dll_node_t* dll_insert(dll_t *list, dll_node_t *cursor, double data);

/*
    Diese Funktion löscht den angegebenen Knoten aus der gegebenen Liste.

    @list (in): Liste, deren Knoten gelöscht werden soll.
    @delnote (in): Knoten, der aus der Liste entfernt werden soll.

    funktioniert soweit.
*/
void dll_delete(dll_t *list, dll_node_t *delnote);

/*
    Diese Funktion fügt einen Wert an den Anfang der Liste ein.

    @list (in): liste, an deren Anfang geschrieben werden soll.
    @data (in): Zahl, die an den Anfagn geschrieben werden soll.

    @list_node (out): Neuer Datenknoten, der nun Am Ende der Liste steht

    soweit getestet.
*/ 
dll_node_t* dll_unshift(dll_t *list, double data);
/*
    Diese Funktion fügt einen neuen Wert ans Ende der Liste ein.

    @list (in): Liste an deren Ende geschrieben werden soll.
    @data (in): Daten, die ans Ende geschrieben werden sollen.

    @list_node (out): Knoten der nun Am Ende der Liste steht.

    soweit getestet.
*/
dll_node_t* dll_push(dll_t *list, double data);

/*
    Diese Funktion entfernt den Wert an erster Stelle der Liste und
    gibt diesen zurück.

    @list (in): liste, deren erstes Element gelöscht werden soll.

    @data  (out): Wert des ersten Elements;

    funktioniert soweit
*/
double dll_shift(dll_t *list);

/*
    Daten am Ende der Liste entfernen und zurückgeben.
    @(dll_t *list)in : double-linked-list pointer

    @data (out): Wert am Ende der Liste.

    fuunktioniert soweit.
*/
double dll_pop(dll_t *list);

/*
    Merge von zwei Listen, wobei die zweite Liste an die erste angefügt wird
    Die beiden Listpointer werden dabei zerstört und nur eine neue gemeinsame
    Liste zurückgegeben.

    @list1 (in): erste Liste, wird zu anfang gestellt.
    @list2 (in): zweite Liste, wird an erste Liste angefügt.

    @out: list-pointer auf gemeinsame Liste.
*/
dll_t* dll_merge(dll_t *list1, dll_t *list2);

/*
    Ausgabe eine Liste,entweder vor- oder rückwärts.

    @list (in): liste die ausgegeben werden  soll.
    @option (in): vor(1)- oder rück(-1)-wärts, Richtung der Ausgabe.

    funktioniert momentan nur vorwärts.
*/
void dll_print(dll_t *list, int option);
#endif