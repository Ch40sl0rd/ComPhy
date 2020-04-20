#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void create_table( double (*fp)(double),
    char *filename,
    double step,
    double start,
    double end
){
    FILE *file =fopen(filename, "w");
    if(!file){
        printf("Fehler bei der Erstellung der  Datei!\n");
        exit(3);
    }
    for(double i = start; i<end; i+=step){
        fprintf(file, "%lf\t%lf\n", i,  fp(i));
    }
    fprintf(file, "%lf\t%lf\n", end, fp(end));
    fclose(file);
}

int main(){
    char* filename  = "cosinus.txt";
    double  start  = 0.;
    double end  = 6.2430;
    double step  = 0.1;
    create_table(cos, filename, step, start, end);
    return  0;
}