#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "arrayhelpers.h"

unsigned long int random(int I_j){
    const int m = 2147483647;
    const int a = 16807;
    const int q = 127773;
    const int r = 2836;
    int I_jp1 = a*(I_j%q) - r*(I_j/q);
    if(I_jp1 <= 0){
        I_jp1 += m;
    }
    return  I_jp1;
}

double rand_frac(int d){
    return (double)d/(double)2147483647;
}

unsigned int index(int i, int j, int length){
    int const i_wrap = (i + length) % length;
    int const j_wrap = (j + length) % length;

    return i_wrap * length + j_wrap;
}

int all_life(unsigned int* array, int length){
    int sum =0;
    for(int i=0; i<length; i++){
        for(int j=0; j<length; j++){
            sum += array[index(i,j,length)];
        }
    }
    return sum;
}

void print_array(unsigned int *array, int l){
    for(int i=0; i<l; i++){
        for(int j=0; j<l-1; j++){
            printf(" %d |", array[index(i, j, l)]);
        }
        printf(" %d\n", array[index(i, l-1, l)]);
    }
}

void swap(unsigned int** a, unsigned int** b){
    unsigned int* temp = *a;
    *a = *b;
    *b = temp;
}

void init_array(unsigned int *array, int l, double threshhold){
    unsigned int state_int = 2345435;
    double state_float = 0.;
    for(int i = 0; i <l; i++){
        for(int j = 0; j<l; j++){
            state_int = random(state_int);
            state_float = rand_frac(state_int);
            if( state_float < threshhold){
                array[index(i, j, l)] = 1;
            }
            else{
                array[index(i, j, l)] = 0;
            }

        }
    }
}

int life_around(unsigned int *array, int length, int x, int y){
    int sum = 0;
    for(int i=-1; i<2;  i++){
        for(int  j=-1; j<2; j++){
            if( i==0 && j==0){
                //do nothing
            }
            else{
                sum+= array[index(x+i, y+j, length)];
            }
        }
    }
    return sum;

}

void iterate_gol(unsigned int *array1, unsigned int *array2,  unsigned int length){
    int life, life_ar, res;
    //Gehe das array1 feld für feld durch
    for(int i=0; i<length; i++){
        for(int j=0; j<length; j++){
            life = array1[i*length +j];
            life_ar = life_around(array1, length, i, j);
            if(life==0){
                if(life_ar==3){
                    res =1;
                }
                else{
                    res=0;
                }
            }
            else if(life==1){
                if((life_ar==2)||(life_ar==3)){
                    res = 1;
                }
                else{
                    res = 0;
                }
            }
            array2[index(i,j,length)] = res;
        }
    }
}

int main(int argc, char *argv[]){
    if(argc < 2){
        printf("Zu wenig optionen. Bitte gebe eine Wahrscheinlichkeit an.\n");
        return -1;
    }
    double const threshhold = atof(argv[1]);

    const unsigned int  l = 12;
    unsigned int Z1[l*l], Z2[l*l];
    unsigned int *Z_tau = Z1;
    unsigned int *Z_taup1 = Z2;

    init_array(Z_tau, l, threshhold);
    print_array(Z_tau, l);
    printf("Die Anzahl der lebenden  Zellen bei %d, %d beträgt %d\n", 11,11, life_around(Z_tau, l, 11,11));
    for(int i=0; i<100; i++){
        iterate_gol(Z_tau, Z_taup1, l);
        swap(&Z_tau, &Z_taup1);
        print_array(Z_tau, l);
        printf("Insgesamt leben noch %d Zellen.\n", all_life(Z_tau,l));
        sleep(1);
        system("clear");
        if(all_life(Z_tau, l) ==  0){
            printf("Alle Zellen sind tot! Es hat %d Spielzüge gebraucht.\n", i+1);
            i=100;
        }
    }
    return 0;
}