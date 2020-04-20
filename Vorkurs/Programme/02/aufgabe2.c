#include <stdio.h>
#include <math.h>

int prim(int c){
    int n = 2;
    while( ((double)n < sqrt((double)c)) && (c%n != 0)){
        n++;
    }
    if(c%n!=0){
        return 1;
    }
    else{
        return 0;
    }
}

int prim_nach(int n){
    int n_0 =  n;
    while(prim(n)==0){
        n++;
    }
    printf("Die nächste Primzahl nach %d ist %d.\n", n_0, n);
    return 0;
}

int main(){
     int list[3] = {20000, 30000, 40000};
     for(int i=0; i<3; i++){
         prim_nach(list[i]);
     }
}