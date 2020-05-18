//make
// ./hausaufgabe2 -n [Anzahl Schritte] -w [Frequenz] -e [obere Grenze]
//Lars Döpper

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#define PI 3.14159265358979323846
#include <math.h>
#include "numerik_own.h"
#include "arrayhelpers.h"
#include "filehelper.h"
#include "h2.h"

void GetUserParam( int argc, char *argv[] );

//declare global variables
double omega_0, gamma, x_0, v_0, end;
int steps;
double alpha, omega_e;

/**
 *  This function represents the vector-valued function of the
 *  differential equation:
 *  x''(t) = -w_0^2*x(t) -2y*x'(t)
 *  This function assumes that the calling methode already
 *  allocated storage for y and f in appropiate dimensions
 * **/
void deq_function(int neq, double t, double *y, double *f){
    if(neq!=2){
        printf("[deq-function] wrong degree of differential equation.\n");
        abort();
    }
    f[0] = y[1];
    f[1] = -omega_0*omega_0*y[0] - 2.0*gamma*y[1];
}

/**
 *  Analytical solution for the differential equation
 *  x''(t) = -w_0^2*x(t) -2y*x'(t)
 *  with starting values x(0) = x_0 and x'(0) = v_0
 *  
 *  \param t: time for the function to be evaluated at.
 * 
 *  \return: value of the function at point t
 * */
double solution(double t){
    double omega_m, pref;
    omega_m = sqrt(-gamma*gamma + omega_0*omega_0);
    pref = (v_0 + x_0*gamma)/omega_m;
    return exp(-gamma*t)*(x_0*cos(omega_m*t) + pref*sin(omega_m*t));
}

/**
 *  Analytical solution for the differential equation
 *  x''(t) = -w_0^2*x(t) -2y*x'(t) + a*cos(w*t)
 *  with starting values x(0) = x_0 and x'(0) = v_0
 *  
 *  \param t: time for the function to be evaluated at.
 * 
 *  \return: value of the function at point t
 * */
double solution2(double t){
    double A, B, w_m, diff_w_sqr, gamma2w_e;
    double pref3, pref4, nenner;
    //Konstanten
    w_m = sqrt(-gamma*gamma + omega_0*omega_0); //Omega aus homogener Lösung
    gamma2w_e = 2.0*gamma*omega_e;  //Faktor aus partikulärer Lösung
    diff_w_sqr = -omega_e*omega_e + omega_0*omega_0;    //Differenz der Quadrate
    //Nenner des Bruches
    nenner = diff_w_sqr*diff_w_sqr + gamma2w_e*gamma2w_e;
    //Vorfaktoren für Funktion.
    A = x_0 - (alpha*diff_w_sqr)/nenner;
    B = v_0/w_m + gamma*A/w_m - (gamma2w_e*omega_e*alpha)/(nenner*w_m);
    pref3 = alpha*diff_w_sqr/nenner;
    pref4 = alpha*gamma2w_e/nenner;
    return exp(-gamma*t)*(A*cos(w_m*t) + B*sin(w_m*t)) + pref3*cos(omega_e*t) + pref4*sin(omega_e*t);
}

/**
 *  This function represents the vector-valued function of the
 *  differential equation:
 *  x''(t) = -w_0^2*x(t) -2y*x'(t) +a*cos(w*t)
 *  This function assumes that the calling methode already
 *  allocated storage for y and f in appropiate dimensions
 * **/
void deq2_function(int neq, double t, double *y, double *f){
    if(neq!=2){
        printf("[deq2-function] wrong degree of differential equation.\n");
        abort();
    }
    f[0] = y[1];
    f[1] = -omega_0*omega_0*y[0] - 2.0*gamma*y[1] + alpha*cos(omega_e*t);
}

/**
 *  This function represents the vector-valued function of the
 *  differential equation:
 *  x''(t) = -w_0^2*x(t) -2y*x'(t) +f_ext(t)
 *  This function assumes that the calling methode already
 *  allocated storage for y and f in appropiate dimensions
 * **/
void deq3_function(int neq, double t, double*y, double*f){
    if(neq!=2){
        printf("[deq2-function] wrong degree of differential equation.\n");
        abort();
    }
    f[0] = y[1];
    f[1] = -omega_0*omega_0*y[0] - 2.0*gamma*y[1] + fext(t);
}

/**
 *  This function  calculates the average power over a number of
 *  of  cycles and returns a list with the average power of each cycle.
 * 
 *  \param omega: frequency of the harmonic oszillator
 *  \param periods: number of periods to analyze.
 * 
 *  \return: list of average powers for each cycle
 * */
double* avg_power(double omega, int periods){
    double h, end, period_duration,sum_power, t, *power;
    double  **dgl_values;
    double *average_power;
    int period_steps;
    //set starting values
    const double start_vals[] = {x_0, v_0};

    //end point for n periods
    end = 2*PI/omega*periods;
    //set omega as new global omega_0
    omega_0 = omega;
    //set step width
    h = 0.001;
    //duration of one period
    period_duration = 2*PI/omega;
    period_steps = (int)(period_duration/h);
    //calculate values for x and x'
    dgl_values = solve_dgl(2, 0.0, end, h, 4, start_vals, &deq3_function);

    average_power = (double*)malloc(sizeof(double)*periods);
    //calculate the average poweer of each cycle
    //iterate over each period
    for(int i=0; i<periods; i++){
        //iterate over every point in a period
        sum_power = 0.0;
        power = (double*)malloc(sizeof(double)*period_steps);
        for(int j=0; j<period_steps; j++){
            t = dgl_values[i*period_steps + j][0];
            //calculate power at every point
            power[j] = fext(t)*dgl_values[i*period_steps+j][2];
            //calculate integrale with trapezoid-rule
        }
        for(int f = 0; f<period_steps-1; f++){
            sum_power += 0.5*h*(power[f] + power[f+1]);
        }
        average_power[i] = sum_power/period_duration;
        free(power);
    }
    
    free2d(dgl_values);
    return average_power;
}

/**
 *  This function calucaltes the average power over a number
 *  of cycles.
 * 
 *  \param omega: frequency to calcualte average for.
 *  \param periods: number of periods to scan for.
 * 
 *  \return average power of all periods.
 * */
double avg_power_sum(double omega, int periods){
    double *powers, sum, avg;
    sum = 0.0;
    avg = 0.0;
    //calculate power for each period seperatly
    powers = avg_power(omega, periods);
    for(int i=0; i<periods; i++){
        //add all powers together
        sum += powers[i];
    }  
    //calculate average power.
    avg = sum/periods;
    free(powers);
    return avg;
}

/**
 *  This function performs a scan with omega and calcutes the power
 *  over 60 cycles for each omega.
 * 
 *  \param start: starting value for omega
 *  \param end: largest value for omega
 *  \param steps: number of steps between start and end
 * 
 * */
void scan_in_omega(double start, double end, int steps){
    double *omegas, *powers;
    omegas = create_array_double(start, end, steps);
    powers = (double*)malloc(sizeof(double)*steps);
    for(int i=0; i<steps; i++){
        powers[i] = avg_power_sum(omegas[i], 60);
    }
    print_data_table_double("Scan_Omega.txt", omegas, powers, steps);
    free(omegas); free(powers);

}

int main( int argc, char *argv[]){
    double h, **val_array, *diff, t;
    double **val_array_2, *diff2;
    double **val_array_3;
    
    GetUserParam(argc, argv);
    //starting values:
    x_0 = 1.5;
    v_0 = 2.75;

    //calculate and initilaze parameters for dgl_solver
    double const start_vals[] = {x_0, v_0};
    h = end/(double)(steps-1);

    //solve first differential equation.
    val_array = solve_dgl(2, 0.0, end, h, 4, start_vals, &deq_function);
    //test rk4 after one step.
    printf("Differnece after one step %15.6e \n", solution(h)-val_array[1][1] );
    printf("Solution at 0: %15.6e\n", solution(0.0));
    /************************************************************************
     * Merke die Index-Notation a[i][j]: i:Zeile, j:Spalte, just for personal use.
     * *********************************************************************/
    printf("Der Wert unserer Funktion bei 0: %15.6e\n", val_array[0][1]);

    diff = create_array_double(0.0, 100.0, steps);
    for(int i=0; i<steps; i++){
        t = val_array[i][0];
        diff[i] = solution(t) - val_array[i][1];
    }
    print_data2file("Runge_Kutta4_diff.txt", diff, steps);

    //Zweite Aufgabe
    //*********************************************************************************************************************************
    alpha = 0.5;
    omega_e = 1.5;

    val_array_2 = solve_dgl(2, 0.0, end, h, 4, start_vals, &deq2_function);
    print_table2file("Kopplung.txt", val_array_2, steps, 3);
    
    diff2 = (double*)malloc(sizeof(double)*steps);
    for(int i=0; i<steps; i++){
        t = val_array_2[i][0];
        diff2[i] = val_array_2[i][1] - solution2(t);
    }
    print_data2file("ExterneKopplung.txt", diff2, steps);
    printf("Exakte Lösung bei 0: %15.6e\n", solution2(0.0));
    //**********************************************************************************************************************************

    //Dritte Aufgabe
    //********************************************************************************************************************************
    //Löse DGL mit neuer Funktion deq3_function
    val_array_3 = solve_dgl(2, 0.0, end, 0.001, 4, start_vals, &deq3_function);
    int steps_3 = end/0.001 + 1;
    print_table2file("ExterneKraft", val_array_3, steps_3, 3);

    double *average_power = avg_power(1.0, 12);
    print_data2file("Power_Average_1.txt", average_power, 12);

    scan_in_omega(1.0, 3.0, 100);


    //*****************************************************************************************************************************
    

    free2d(val_array);  free2d(val_array_2); free2d(val_array_3);
    free(diff2);
    free(diff);
    free(average_power);
    return 0;
}


void GetUserParam( int argc, char *argv[] ){
/* Variablen: */
    int i;
    char* endptr;
    const char usage[] = 
        "# hausaufgabe2 [-n <Schrittzahl>] [-w <Frequenz>] [-e <Obere Grenze>]";
    const char error_message[] =
        "# FEHLER(GetuserParam): falsche Option: ";

    if (argc>1) { /* falls es ueberhaupt Parameter gibt ... */
        for (i=1; i<argc; i++){
            /* parameter 2 Charakter lang und sollte mit '-' anfaengen ... */
            if ( (strlen(argv[i])==2) && (argv[i][0] == '-') ) { 
            switch (argv[i][1]) { 
		    case 'w':
                omega_0 = strtod( argv[++i], &endptr);
                if ( (!isspace(*endptr) && (*endptr) != 0) ) {
			        printf(" %s \n %s \n",error_message,usage);
                    exit(1);
                }
                break;
		    case 'e':
                end = strtod( argv[++i], &endptr);
                if ( (!isspace(*endptr) && (*endptr) != 0) ) {
			        printf(" %s \n %s \n",error_message,usage);
			        exit(1);
                    }
                break;
            case 'n':
                steps = strtol( argv[++i], &endptr, 10);
                if ( (!isspace(*endptr) && (*endptr) != 0) ) {
			        printf(" %s \n %s \n",error_message,usage);
			        exit(1);
                }
                break;
            default:
		        printf(" %s \n %s \n",error_message,usage);
		            exit(1);
                }
            } 
            else {
	            printf(" %s \n %s \n",error_message,usage);
	            exit(1);
            } /* end of: if-then-else */
        } /* end-of: for */ 
        gamma = 0.1;
    } /* end-of: if */
    else
    {      
        printf("[GetUserParams] No values declared, use default values.\n");
        //set default values
        omega_0 = 1.0;
        gamma = 0.1;
        steps = 1000;
        end = 5.0;
    }  
}