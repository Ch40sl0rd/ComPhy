#include <stdlib.h>
#include "arrayhelpers.h"

int* merge(int *list1, int n1, int *list2, int n2){
    int pos1, pos2, *newlist, i;
    newlist = (int*)malloc(sizeof(int)*(n1+n2));
    if(!newlist){
        printf("ERROR bei Speicherallokierung!\n");
        exit(1);
    }
    //Gehe die beiden Liste elementeweise durch und füge das kleinste ein.
    for(i=0, pos1=0, pos2=0; (pos1<n1)&&(pos2<n2) ; i++){
        if(list1[pos1] <= list2[pos2]){
            newlist[i] = list1[pos1];
            pos1++;
        }
        else if(list1[pos1] > list2[pos2]){
            newlist[i] = list2[pos2];
            pos2++;
        }
        else{
            printf("Das sollte nicht passieren!\n");
        }
    //Der Rest wird überschrieben  
    }
    if( (pos1==n1) && (pos2<n2) ){
        for(;pos2<n2; i++){
            newlist[i] = list2[pos2];
            pos2++;
        }
    }
    if( (pos2==n2)&&(pos1<n1)){
        for(;pos1<n1; i++){
            newlist[i] = list1[pos1];
            pos1++;
        }
    }
    return newlist;
}

void mergesort(int *list, int n){
    int *newlist;
    int border;
    //Wir erstellen eine Grenze zwischen den zwei Listen, je nach länge
    if(n<=1){
        return;
    }
    if(n%2==0){
        border = n/2;
    }
    else{
        border = (n+1)/2;
    }

    //Mergesort auf beide Unterlisten ausführen
    mergesort(list, border);
    mergesort(&(list[border]), n-border);

    //Eine Neue Liste aus den beiden sortiereten Listen erstellen.
    newlist = merge(list, border, &(list[border]), n-border);
    for(int i=0; i<n; i++){
        list[i] = newlist[i];
    }
    free(newlist);


}

int main(int argc, char* argv[]){
    if(argc<2){
        printf("Bitte gebe eine Datei an, die du sortieren möchtest.\n");
        exit(2);
    }

    int *list, n, dat;

    const char* dateiname = argv[1];
    FILE  *file = fopen(dateiname, "r+");
    //Gehe die Datei Zeile für Zeile durch und merke dir die Zeilennummer
    int i=0;
    //feof() gibt zurück, ob der Dateipointer schon am Ende angekommen ist.
    while(!feof(file)){
        fscanf(file, "%d\n", &dat);
        //printf("Test  %d\n", i);
        i++;
    }
    n=i;
    rewind(file);

    //Speichere  die Liste lokal ab.
    list = (int*)malloc(sizeof(int)*n);
    if(!list){
        printf("Error bei Speicherallokierung!\n");
        exit(1);
    }

    for(i=0; i<n;  i++){
        fscanf(file, "%d\n", &dat);
        list[i] = dat;
    }
    
    array_ausgabe(list, n, 2);
    mergesort(list, n);
    array_ausgabe(list, n, 2);
    rewind(file);
    for(i=0; i<n;  i++){
        fprintf(file, "%d\n", list[i]);
    }

    //Freigabe des Speciheers und des Dateizeigers.
    free(list);
    fclose(file);
    return 0;
}