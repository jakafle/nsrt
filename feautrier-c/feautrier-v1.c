/* Feautrier Method for radiation transfer */

/* Use the standard and GSL libraries */
#include <stdio.h>
#include <gsl/gsl_linalg.h>

/* Declaration of source function 
Takes two arguments: optical depth (tau) and flux value (flux_val) */
double source_func(double tau, double flux_val);

/* The main function, declare array in a heap and update it using a for loop */
int main(void){
    
    const int N = 6; // take six layers of optical depth
    double mu = 0.5773502691896257; // gaussian quadrature for angles
    double deltau = 0.1; // the step size for optical depth

    /* create vectors for optical depth, feautrier variable u, and source values */
    double tau[N];
    double svals[N];
    double evals[N];
 
    /* Implement the upper boundary conditions */
    double a = pow(mu,2)/pow(deltau,2);
    double d = 1 + (pow(mu,2)*2)/(pow(deltau,2));
    double c = pow(mu,2)/pow(deltau,2);

    /* Implement the lower boundary conditions */
    double d0 = 1 + mu/deltau + deltau/(2*mu);
    double c0 = mu/deltau;
    double d5 = 0.5 + (mu/deltau);
    double a5 = (mu/deltau) - 0.5;
    
    /* Create the matrix */
    double amat[] = {d0, -c0, 0, 0, 0, 0, -a, d, -c, 0, 0, 0, 0, -a, d, -c, 0, 0, 0, 0, -a, d, -c, 0, 0, 0, 0, -a, d, -c, 0, 0, 0, 0, -a5, d5};

    /* Update the value of optical depth */
    for (int i=1; i < N; i++)
    {
        tau[i] = tau[i-1] + 0.1;
    }

    /* Find the evals from the source function */
    for (int i=0; i < N; i++)
    {
        evals[i] = source_func(tau[i], 1);
    }

    /* Find the right side of the matrix equation by implementing the boundary conditions */    
    svals[0] = (evals[0]*deltau)/(2*mu);
    svals[1] = evals[1];
    svals[2] = evals[2];
    svals[3] = evals[3];
    svals[4] = evals[4];
    svals[5] = (evals[5]*(0.5 + (mu/deltau)) + evals[4]*(0.5 - (mu/deltau)));

    /* Using GSL solving a matrix equation 
        Use LU decomposition to solve Au = b*/
    gsl_matrix_view A = gsl_matrix_view_array(amat, 6, 6);
    gsl_vector_view b = gsl_vector_view_array (svals, 6);

    gsl_vector *u = gsl_vector_alloc (6);

    int s;

    gsl_permutation * p = gsl_permutation_alloc (6);

    gsl_linalg_LU_decomp (&A.matrix, p, &s);

    gsl_linalg_LU_solve (&A.matrix, p, &b.vector, u);

    printf ("The feautrier variable is u = \n");
    gsl_vector_fprintf (stdout, u, "%g");

    gsl_permutation_free (p);
    gsl_vector_free (u);
    
    return 0;
}

/* The source function */
double source_func(double tau, double flux_val)
{
    double bval = (3./4)*flux_val*(tau + 2./3);
    return bval;
}
