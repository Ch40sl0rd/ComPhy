#include <stdio.h>
#include <math.h>

int ggT(int a, int b){
    while(1){
        if(a==0){
            return b;
        }
        if(b==0){
            return a;
        }
        if(a>b){
            a = a%b;
        }
        else{
            b = b%a;
        }
    }

}
int main(){
    int m=5498384,  n=215844;
    printf("Der GgT von %d und %d ist %d.\n", m,n, ggT(m,n));
    return 0;
}