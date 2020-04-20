#include<stdlib.h>
#include<stdio.h>
int main(){
   int n,i;
   double *angle;
   double *cosangle;
   FILE *input=fopen("input","r");
   if (input == NULL){
     printf("Fehler, wir konnen leider nicht die input datei eröffnen\n");
     exit(1);
   } 
   //Wir lesen die Anzahl der folgende Zeilen Aus
   fscanf(input,"%d\n", &n);
   printf("Wir haben %d Zeilen im input Datei\n",n);
   angle=(double*)malloc(sizeof(double)*n);
   cosangle=(double*)malloc(sizeof(double)*n);
   if ( (angle == NULL) && (cosangle == NULL)){
     printf("Wir haben nicht genug Speichern");
     exit(1);
   }
   for (i=0;i<n;++i){
     fscanf(input,"cos(%le) = %le\n",&angle[i],&cosangle[i]);
   }
   for (i=0;i<n;++i)
     printf("cos(%le)=%le\n",angle[i],cosangle[i]);
   free(angle);
   free(cosangle);
   fclose(input);
}
