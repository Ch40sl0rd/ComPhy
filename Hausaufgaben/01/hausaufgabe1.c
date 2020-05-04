// make
// ./hausaufgabe1
// Lars Döpper
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define M_PI acos(-1.)
#include "arrayhelpers.h"
#include "filehelper.h"
#include "numerik_own.h"

/**
 *  This function calculates the value of the integrand
 *  of the electric potential. You can use -1000 and 1000 as borders,
 *  the value of these borders equates to 0 in double precission.
 * 
 *  \param x: x value for integration
 *  \param p: list of all paramerts
 * 
 *  \return: value of the integrand at x.
 * */
double integrand(double x, void* p){
    double a,z;
    a = (double)((double*)p)[0];
    z = (double)((double*)p)[1];
    return exp(-(x*x)/(a*a))/(sqrt(x*x + z*z));
}

double integrand_2(double x, void*p){
    double a,z;
    a = (double)((double*)p)[0];
    z = (double)((double*)p)[1];
    return exp(-(x*x)/(a*a))/pow(x*x + z*z, 1.5);
}

double electric_field(double z, void*p){
    double a, pre_f, Q, alpha; //constant factors
    double integrale; //value of the integrale
    double res; //result
    double acc;
    //reading all constant factors from list of parameters.
    a = (double)((double*)p)[0];
    pre_f = (double)((double*)p)[1];
    Q = (double)((double*)p)[2];
    alpha = (double)((double*)p)[3];
    acc = (double)((double*)p)[4];

    //creating list of parameters for integrand function
    double *p_in = (double*)malloc(sizeof(double)*2.);
    if(!p){
        printf("Error whiel allocating memory.\n");
        exit(1);
    }
    p_in[0] = a;
    p_in[1] = z;
    integrale = integrate_romberg(-20., 20., 100, &integrand_2, acc, 15, (void*)p_in);
    res = alpha*Q*pre_f/(a)*z*integrale;
    free(p_in);
    return res;
}

/**
 *  This function calculates the value of the electric field of two
 *  charges with distance d.
 * 
 *  \param z: point to evaluate the electric field at.
 *  \param p: list of all parameters.
 * 
 *  \return: vlaue of the elctric field at z.
 * */
double electric_field_dual(double z, void* p){
    //declare all parameters:
    double a1, Q1, a2, Q2, d;
    double alpha, pre_f, acc;
    double res;
    double p_in[5];
    //read values from list of parameters
    a1 = (double)((double*)p)[0];
    a2 = (double)((double*)p)[1];
    Q1 = (double)((double*)p)[2];
    Q2 = (double)((double*)p)[3];
    pre_f = (double)((double*)p)[4];
    alpha = (double)((double*)p)[5];
    acc = (double)((double*)p)[6];
    d = (double)((double*)p)[7];

    //create list of parameters for electric field function:
    p_in[0] = a1;
    p_in[1] = pre_f;
    p_in[2] = Q1;
    p_in[3] = alpha;
    p_in[4] = acc;
    res = electric_field(z, (void*)p_in);

    //swap parameters fro second part
    p_in[0] = a2;
    p_in[2] = Q2;
    res -= electric_field((d-z), (void*)p_in);

    return res;
}

/**
 *  This function calculates the given electrostatic potential by
 *  integrating the given function
 *  \param z: position of the test charge
 *  \param p: pointer with all the parameters for the function.
 *  \param acc: accuracy of the integration.
 * 
 *  \return: electrostatic potential at point z.
 **/
double potential(double z, void* p){
    double a, pre_f, Q, alpha; //constant factors
    double integrale; //value of the integrale
    double res; //result
    double acc;
    //reading all constant factors from list of parameters.
    a = (double)((double*)p)[0];
    pre_f = (double)((double*)p)[1];
    Q = (double)((double*)p)[2];
    alpha = (double)((double*)p)[3];
    acc = (double)((double*)p)[4];

    //creating list of parameters for integrand function
    double *p_in = (double*)malloc(sizeof(double)*2.);
    if(!p){
        printf("Error whiel allocating memory.\n");
        exit(1);
    }
    p_in[0] = a;
    p_in[1] = z;
    integrale = integrate_romberg(-20., 20., 100, integrand, acc, 15, (void*)p_in);
    //printf("Wert des Integrals: %15.6e\n", integrale);
    res = pre_f * (alpha*Q/a)*integrale;
    //printf("Wert des Ergebnisses: %15.6e\n", res);
    free(p_in);
    return res;
}

/**
 *  
 * 
 * 
 * 
 * */
double potential_dual(double z, void* p){
    //declare all constant paramerters.
    double a1, a2, Q1, Q2, pre_f, alpha, d;
    double potential1, potential2;
    double res;
    double acc;
    double p_in[5];
    //read all data from the list of parameters:
    a1 = (double)((double*)p)[0];
    a2 = (double)((double*)p)[1];
    Q1 = (double)((double*)p)[2];
    Q2 = (double)((double*)p)[3];
    pre_f = (double)((double*)p)[4];
    alpha = (double)((double*)p)[5];
    acc = (double)((double*)p)[6];
    d = (double)((double*)p)[7];
    //create list of parameters for fist potential
    p_in[0] = a1;
    p_in[1] = pre_f;
    p_in[2] = Q1;
    p_in[3] = alpha;
    p_in[4] = acc;
    potential1 = potential(z, (void*)p_in);
    //switch parameters for second potential
    p_in[0] = a2;
    p_in[2] = Q2;
    potential2 = potential((d-z), (void*)p_in);

    res = potential1 + potential2;
    return res;
}

int main(){
    //List of all paramerts needed for this homework
    double Q_1, Q_2, a_1[6], p2[8], pre_f; //all constant coeffiecents
    double acc; //accuracy of the integration.
    double * z, *z2, *pot_1, *e_field,*e_field2, *pot_2, *e_field_calc, *e_field2_calc, *a_2, *zero_points;
    int length; //number od datapoints.
    double a1, a2, b1, b2; //borders of the potential.

    a_1[0] = 4.;
    a_1[1] = 2.;
    a_1[2] = 1.;
    a_1[3] = 0.1;
    a_1[4] = 0.01;
    a_1[5] = 0.001;
    

    char aufgabe1[6][20] = {"a_4.txt", "a_2.txt", "a_1.txt", "a_0_1.txt", "a_0_01.txt", "a_0_001.txt"};
    //Coeffitients for inetgration
    length = 1000;
    a1 = 0.;
    b1=5.;

    Q_1 = 1.;
    pre_f = 1.;
    acc = 1e-8;

    //create list of parametrs for potential.
    double *p_1 = (double*)malloc(sizeof(double)*5);
    if(!p_1){
        printf("Error while allocating memory.\n");
        exit(1);
    }
    p_1[1] = pre_f;
    p_1[2] = Q_1;
    p_1[3] = 1/sqrt(M_PI);
    p_1[4] = acc;

    z = create_array_double(a1, b1, length);
    //calculate potential for diffenrent values of a_1
    printf("Wir erstellen jetzt die Potentiale für verschiedene Werte von a_1:\n\n");
    for(int i=0; i<6; i++){
        p_1[0] = a_1[i];
        pot_1 = create_array_function_double(z, potential, length, (void*)p_1);
        print_data_table_double( aufgabe1[i] , z, pot_1, length);
        free(pot_1);
    }

    //creating the electric field with a=1.
    p_1[0] = a_1[2];
    pot_1 = create_array_function_double(z, potential, length, (void*)p_1); //calculating potential
    //use the first dervivative for the elekktric field
    printf("Wir berechnen jetzt das Elektrische Feld der Ladungsverteilung aus numerischer Ableitung.\n\n");
    e_field = derivate_sym_one_array(pot_1, a1, b1, length);
    //richtiger Vorfaktor
    for(int i=0; i<length;i++){
        e_field[i] *=-1.;
    } 
    print_data_table_double( "E-Feld.txt" , z, e_field, length);
    printf("Wir berechnen jetzt das Elektrische Feld der Ladungsverteilung mit numerischer Integration:\n\n");
    e_field_calc = create_array_function_double(z, electric_field, length, (void*)p_1);
    print_data_table_double("E_Feld_calc.txt", z, e_field_calc, length);


    //calculating the double potential
    Q_2 = 4.;
    a2 = 0.;
    b2 = 12.;
    z2 = create_array_double(a2, b2, length);
    //create list of parameters:

    p2[0] = 4.;
    p2[1] = 4.;
    p2[2] = Q_1;
    p2[3] = Q_2;
    p2[4] = pre_f;
    p2[5] = 1/sqrt(M_PI);
    p2[6] = acc;
    p2[7] = b2;
    printf("Wir erstellen jetzt das Potential der doppelten Ladungsverteilung:\n\n");
    pot_2 = create_array_function_double(z2, &potential_dual, length, p2);
    print_data_table_double("Pot_Dual.txt", z2, pot_2, length);
    printf("Wir berechnen jetzt numerisch das E-Feld der doppelten Ladungsverteilung:\n\n");
    e_field2 = derivate_sym_one_array(pot_2, a2, b2, length);
    for(int i=0; i<length; i++){
        e_field2[i] *= -1.;
    }
    print_data_table_double("e_field_2.txt", z2, e_field2, length);
    e_field2_calc = create_array_function_double(z2, electric_field_dual, length, (void*)p2);
    print_data_table_double("e_field2_calc.txt", z2, e_field2_calc, length);

    printf("Die Nullstelle des E-Feldes lautet: %15.6e\n", zero_crossing(&electric_field_dual, 11.5, 11., acc, (void*)p2));

    //erstelle Werte für die Nullstellenbestimmung des elektrischen Feldes.
    a_2 = create_array_double(0.4, 4., 100);
    zero_points = (double*)malloc(sizeof(double)*100);
    if(!zero_points){
        printf("Error while allocating memory.\n");
        exit(1);
    }
    for(int i=0; i<100; i++){
        p2[1] = a_2[i];
        zero_points[i] = zero_crossing(&electric_field_dual, 11.5, 11., acc, (void*)p2);
    }

    print_data_table_double("Nullstellen.txt", a_2, zero_points, 100);


    //free all used memory.
    free(z); free(pot_1); free(p_1); free(e_field);free(e_field_calc); free(pot_2);
    free(z2); free(e_field2); free(e_field2_calc); free(a_2); free(zero_points);
}