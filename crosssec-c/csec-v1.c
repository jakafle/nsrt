#include <stdio.h> 
#include <math.h> 
#include <stdlib.h>
#include <complex.h>

/*---------------------- Constants used throughout the codes --------------------*/

#define CHARGE_E 1.602176634e-19        // electron charge in Coulomb
#define PLANCK_CONST 6.62607015e-34     // Planck constant (Js)
#define SPEED_LIGHT 299792458.0         // Speed of light (m/s)
#define EPSILON_0 8.8541878128e-12      // vacuum permittivity (epsilon_0)
#define MASS_E 9.1093837015e-31         // mass of electron (kg)
#define PI 3.141592653589793            // pi
#define RADIUS_E 2.8179403262049284e-15 // classical electron radius (r_e) in m

/*------------------------ Function Prototypes ----------------------------------*/

double energy_to_freq(double E); // E is in keV
double complex b_calc(double theta); // theta is in radians


/*------------------------ complex datatype -------------------------------------*/
typedef struct complexn
{
    double real;
    double imag;
} complexn;

/*----------------------- Variables throughout the code -------------------------*/
// define density?

int main(void)
{   
    // initial frequency
    double E_i = 30.0; // in keV
    double omega_i = energy_to_freq(E_i); // frequency for initial energy

    // radiation damping
    double rad_damp = (2*pow(CHARGE_E,2)*pow(omega_i,2))/(2*MASS_E*pow(SPEED_LIGHT,3));

    // cyclotron frequency
    double E_c = 20.0; // in keV
    double omega_c = energy_to_freq(E_c);

    // wave number
    int N = 20; // b parameter array size

    // allocate memory for the array
    complexn *b_param = (complexn *)malloc(N * sizeof(complexn));

    // check memory allocation
    if (b_param == NULL) return 1;

    for (int i = 0; i < 90; i + 10)
    {
        b_param[i].real = creal(b_calc(i));
        b_param[i].imag = cimag(b_calc(i));
    }

    return 0;
}

// function to convert energy to angular frequency
double energy_to_freq(double E)
{
    // E = hf => f = E/h
    double omega = 2*PI*(E*1000*CHARGE_E)/PLANCK_CONST;
    return omega;
}

// function to return the b parameter
double complex b_calc(double theta)
{
    // convert theta to radians
    double angle_rad = theta*(PI/180);
    return 0.0;
}
