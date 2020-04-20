#include "bucketsort.h"

int bucket_sort(unsigned short *array, unsigned len) {

	unsigned *buckets, max, i, array_index;
	unsigned short k;

	/* Finde das Maximum des arrays. */
	for (i=0; i<len; i++)
		if (array[i] > max) max = array[i];

	/* Da wir auch in array[max] mindestens einmal
	 * etwas eintragen wollen, müssen wir eine
	 * Speicherzelle mehr als max allokieren. */
	max += 1;

	/* Allokiere Speicher für die Buckets. */
	if ( (buckets = malloc(max*sizeof(int))) ) {
		/* Alle Buckets zu 0 initialisieren. */
		memset(buckets, 0, max*sizeof(int));
		/* Für jeden eintrag im array erhöhen wir nun den
		 * Wert des entsprechenden Buckets um eins. */
		for (i=0; i<len; i++)
			buckets[array[i]]++;
		/* Wir müssen nun den Array neu beschreiben: */
		array_index = 0;
		/* Für jede mögliche Zahl, die vorkommen kann, */
		for (k=0; k<max; k++)
			/* iteriere durch den entsprechenden Bucket
			 * und fülle array langsam auf. */
			for (i=0; i<buckets[k]; i++)
				array[array_index++] = k;
		/* Wir können nun den Speicher für die Buckets
		 * wieder freigeben. */
		free(buckets);
		return 0;
	} else {
		return 1;
	}
}

