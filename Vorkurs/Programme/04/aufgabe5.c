#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include "arrayhelpers.h"
#include "mathe_own.h"

void sort_lang(int*  array, int length){
    int small;
    int pos;
    for(int i=0; i<length; i++){
        small = array[i];
        pos  = i;
        for(int j=i+1; j<length; j++){
            if(array[j] <= small){
                small = array[j];
                pos  = j;
            }
        }
        swap_pos(array, length, i, pos);
    }
}

void bucketsort(int* array,  int length){
    int max  = array[0];
    
    for(int i = 0; i<length; i++){
        if(max < array[i]){
            max = array[i];
        }
    }
    max+=1;
    //Wir benötigen eine leere Liste mit allen möglichen Einträgen.
    int* buckets;
    int array_index;
    if( (buckets = malloc(sizeof(int)*max)) ){
        memset(buckets, 0, max*sizeof(int));
        for(int i = 0; i<length; i++){
            buckets[array[i]]++;
        }
        //Jetzt überschreiben wir das array mit den Werten aus den Buckets.
        //arrayindex läuft hier von 0 bis length durch um jedes element von array ansprechen zu können.
        array_index = 0;
		/* F�r jede m�gliche Zahl, die vorkommen kann, */
		for (int k=0; k<max; k++)
			/* iteriere durch den entsprechenden Bucket
			 * und f�lle array langsam auf. */
			for (int i=0; i<buckets[k]; i++)
				array[array_index++] = k;
		/* Wir k�nnen nun den Speicher f�r die Buckets
		 * wieder freigeben. */
		free(buckets);
    }
    else{
        printf("Error bei malloc.\n");
    }

}

int main(){
    const int n = 50;
    int array[n];
    array_init_std(array, n);
    array_ausgabe(array, n, 2);
    rotate_k(array, n, 24);
    array_ausgabe(array, n, 2);
    bucketsort(array, n);
    array_ausgabe(array, n, 2);
}