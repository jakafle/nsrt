// include the libraries
#include <stdio.h>
#include <math.h>

#define KB 1.38e-23
#define PI 3.141592654
#define HC 1.05e-34
#define R 2.81e-15

/*
The function takes initial frequency w, final frequency w_prime, initial angle theta,
and final angle theta_prime

w_prime will be a range of frequency values that the redistribution can happen,
theta_prime will be chosen to be a few angles where the photon can direct to after scattering
*/
double crosssec(double w, double w_prime, double theta, double theta_prime);

// the maxwell-boltzmann distribution
double maxwelldist(double p, double T);

// pi_plus
double pi_plus(double wCyc, double w, double w_prime, double theta, double theta_prime);

// pi_minus
double pi_minus(double wCyc, double w, double w_prime, double theta, double theta_prime); 

// pi_z
double pi_z(double wCyc, double w, double w_prime, double theta, double theta_prime);

// Function to perform numerical integration using Gaussian Quadrature
double gaussian_quadrature(double (*func)(double, double, double, double), double lower, double upper, int n, double E, double E_prime, double me) {
  // Pre-computed weights and abscissas for Gaussian Quadrature (replace with your table)
  double *weights = malloc(n * sizeof(double)); 
  double *abscissas = malloc(n * sizeof(double));

  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    // Transform abscissas to the integration interval
    double x = (lower + upper) / 2.0 + (upper - lower) / 2.0 * abscissas[i];
    sum += weights[i] * func(x, E, E_prime, me);
  }

  free(weights);
  free(abscissas);
  return (upper - lower) / 2.0 * sum;
}

int main() {

  // Define integrand function (replace with your actual f(p) function)
  double (*integrand_ptr)(double, double, double, double) = maxwelldist;

  // print the cross-section values in a loop for different frequencies and angles
  // cross_sec = pow(R,2)*(wprime/w)*(m/del_k)*maxwelldist*(pi_plus + pi_minus + pi_z);
  return 0;
}

double maxwelldist(double p, double T) {
    // return the maxwellian
    double f_val = pow((2*PI*KB*T),-0.5)*exp(-pow(p,2)/(2*KB*T));

    return f_val;
}

// pi_plus
double pi_plus(double wCyc, double w, double w_prime, double theta, double theta_prime) {

    double kprime = w_prime*cos(theta_prime);
    double k = w*cos(theta);
    double del_k = kprime - k;
    double del_w = w_prime - w;
    double p_c = (del_w/del_k) + HC*del_k*0.5;
    double pi_p = 1 - (wCyc/(w_prime + wCyc - p_c*kprime + 0.5*pow(kprime,2)));

    return pi_p;
}

// pi_minus
double pi_plus(double wCyc, double w, double w_prime, double theta, double theta_prime) {

    double kprime = w_prime*cos(theta_prime);
    double k = w*cos(theta);
    double del_k = kprime - k;
    double del_w = w_prime - w;
    double p_c = (del_w/del_k) + HC*del_k*0.5;
    double pi_m = 1 - (wCyc/(w_prime - wCyc - p_c*k - 0.5*pow(k,2)));

    return pi_m;
}

// pi_z
double pi_z(double wCyc, double w, double w_prime, double theta, double theta_prime) {

    double kprime = w_prime*cos(theta_prime);
    double k = w*cos(theta);
    double del_k = kprime - k;
    double del_w = w_prime - w;
    double p_c = (del_w/del_k) + HC*del_k*0.5;
    double pi_z = 1 + ((p_c + 0.5*k)*(p_c + k - 0.5*kprime))/(w - p_c*k - 0.5*pow(k,2)) - ((p_c - 0.5*kprime)*(p_c - kprime + 0.5*k))/(w_prime - p_c*kprime - 0.5*pow(kprime,2));

    return pi_z;
}