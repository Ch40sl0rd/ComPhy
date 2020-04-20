#include <stdio.h>
#include <math.h>

double power_own(double x, int n){
    if(n==0){
        return 1;
    }
    if(n%2==0){
        return power_own(x*x, n/2);
    }
    else{
        return x*power_own(x*x, (n-1)/2);
    }
}

double power_for(double x, int n){
    //eine naive implementierung einer potenzfunktion.
    double ergebnis = 1.;
    for(int i=0; i<n;i++){
        ergebnis*=x;
    }
    return ergebnis;
}

double power_own_for(double x, int n){
    /*
    Wir benötigen hier  einen neuen Zwischenspeicher xi für die momentane Basis. Sonst verändert sich nicht viel von der Methode oben.
    */
    double ergebnis = 1.;
    double xi = x;
    int i = 0;
    for(i=n; i>0;){
        if(i%2==0){
            i = i/2;
            xi = xi*xi;
        }
        else if(i%2!=0){
            i = (i-1)/2;
            ergebnis = xi*ergebnis;
            xi = xi* xi;
        }
    }
    return ergebnis;
}

int main(){
    double x = 0.9999999999;
    int n = 2000000000;
    double x_n = 3.;
    int n_n = 3;
    printf("Teste die Exponentialfunktion: %f\n", power_own(x,n));
    printf("Teste die neue Exponentialfunktion: %f\n", power_for(x,n));
    printf("Teste die Exponentialfunktion: %f\n", power_for(x_n,n_n));
    printf("Teste die Exponentialfunktion: %f\n", power_own_for(x,n));
    return 0;
}