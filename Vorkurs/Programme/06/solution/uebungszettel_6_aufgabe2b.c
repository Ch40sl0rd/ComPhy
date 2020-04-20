#include<stdlib.h>
#include<stdio.h>
int main(){
   int n,i;
   double angletmp;
   double cosangletmp;
   double *angle;
   double *cosangle;
   FILE *input=fopen("input","r");
   if (input == NULL){
     printf("Fehler, wir konnen leider nicht die input datei eröffnen\n");
     exit(1);
   } 
   //Wir lesen die Anzahl der folgenden aus
   for (i=0;fscanf(input,"cos(%le) = %le\n",&angletmp,&cosangletmp)!=EOF;++i);
   printf("Anzahl der Eintraege %d\n", i);
   n=i;
   rewind(input);
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
