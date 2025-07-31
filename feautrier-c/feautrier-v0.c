#include <stdio.h>
#include <gsl/gsl_linalg.h>

double source_func(double tau, double flux_val);

int main(void){
    
    const int N = 6;
    double mu = 0.5773502691896257;
    double deltau = 0.1;
    double tau[N];
    double svals[N];
    double evals[N];
    
    double a = pow(mu,2)/pow(deltau,2);
    double d = 1 + (pow(mu,2)*2)/(pow(deltau,2));
    double c = pow(mu,2)/pow(deltau,2);

    double d0 = 1 + mu/deltau + deltau/(2*mu);
    double c0 = mu/deltau;
    double d5 = 0.5 + (mu/deltau);
    double a5 = (mu/deltau) - 0.5;
    
    double amat[] = {d0, -c0, 0, 0, 0, 0, -a, d, -c, 0, 0, 0, 0, -a, d, -c, 0, 0, 0, 0, -a, d, -c, 0, 0, 0, 0, -a, d, -c, 0, 0, 0, 0, -a5, d5};
        
    for (int i=0; i <= N; i++)
        {
            tau[i] = tau[i-1] + 0.1;
        }

    for (int i=0; i <= N; i++)
        {
            evals[i] = source_func(tau[i], 1);
        }

    svals[0] = (evals[0]*deltau)/(2*mu);
    svals[-1] = (evals[-1]*(0.5 + (mu/deltau)) + evals[-2]*(0.5 - (mu/deltau)));
    svals[1] = evals[1];
    svals[2] = evals[2];
    svals[3] = evals[3];
    svals[4] = evals[4];

    gsl_matrix_view m = gsl_matrix_view_array(amat, 6, 6);
    gsl_vector_view b = gsl_vector_view_array (svals, 6);

    gsl_vector *x = gsl_vector_alloc (6);

    int s;

    gsl_permutation * p = gsl_permutation_alloc (6);

    gsl_linalg_LU_decomp (&m.matrix, p, &s);

    gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

    printf ("x = \n");
    gsl_vector_fprintf (stdout, x, "%g");

    gsl_permutation_free (p);
    gsl_vector_free (x);
    
    return 0;
}

double source_func(double tau, double flux_val)
{
    double bval = (3/4)*flux_val*(tau + 2/3);
    return bval;
}