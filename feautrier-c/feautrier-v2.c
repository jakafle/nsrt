#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>

/* Declaration of source function 
Takes two arguments: optical depth (tau) and flux value (flux_val) */
double source_func(double tau, double flux_val);

int main(void) {
    
    int N = 6; // Number of layers of optical depth
    double mu = 0.5773502691896257; // Gaussian quadrature for angles
    double deltau = 0.1; // Step size for optical depth

    /* Create vectors for optical depth, Feautrier variable u, and source values */
    double *tau = malloc(N * sizeof(double));
    double *svals = malloc(N * sizeof(double));
    double *evals = malloc(N * sizeof(double));

    /* Implement the upper boundary conditions */
    double a = mu * mu / (deltau * deltau);
    double d = 1 + 2 * mu * mu / (deltau * deltau);
    double c = mu * mu / (deltau * deltau);

    /* Implement the lower boundary conditions */
    double d0 = 1 + mu / deltau + deltau / (2 * mu);
    double c0 = mu / deltau;
    double dN_minus_1 = 0.5 + (mu / deltau);
    double aN_minus_1 = (mu / deltau) - 0.5;
    
    /* Create the matrix */
    double *amat = malloc(N * N * sizeof(double));
    
    /* Update the value of optical depth */
    for (int i = 0; i < N; i++) {
        tau[i] = i * deltau;
    }

    /* Find the evals from the source function */
    for (int i = 0; i < N; i++) {
        evals[i] = source_func(tau[i], 1);
    }

    /* Find the right side of the matrix equation by implementing the boundary conditions */    
    svals[0] = (evals[0] * deltau) / (2 * mu);
    svals[N - 1] = (evals[N - 1] * dN_minus_1 + evals[N - 2] * (0.5 - (mu / deltau))) / (2 * mu);

    /* Construct the tridiagonal matrix dynamically */
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                if (i == 0) {
                    amat[i * N + j] = d0;
                } else if (i == N - 1) {
                    amat[i * N + j] = dN_minus_1;
                } else {
                    amat[i * N + j] = d;
                }
            } else if (i == j - 1) {
                if (i == N - 2) {
                    amat[i * N + j] = -aN_minus_1;
                } else {
                    amat[i * N + j] = -c;
                }
            } else if (i == j + 1) {
                if (i == 0) {
                    amat[i * N + j] = -c0;
                } else {
                    amat[i * N + j] = -a;
                }
            } else {
                amat[i * N + j] = 0.0;
            }
        }
    }

    /* Using GSL to solve a matrix equation (LU decomposition) */
    gsl_matrix_view A = gsl_matrix_view_array(amat, N, N);
    gsl_vector_view b = gsl_vector_view_array(svals, N);

    gsl_vector *u = gsl_vector_alloc(N);

    int s;
    gsl_permutation *p = gsl_permutation_alloc(N);

    gsl_linalg_LU_decomp(&A.matrix, p, &s);
    gsl_linalg_LU_solve(&A.matrix, p, &b.vector, u);

    printf("The Feautrier variable is u = \n");
    gsl_vector_fprintf(stdout, u, "%g");

    /* Free allocated memory */
    free(tau);
    free(svals);
    free(evals);
    free(amat);
    gsl_permutation_free(p);
    gsl_vector_free(u);

    return 0;
}

/* The source function */
double source_func(double tau, double flux_val) {
    double bval = (3. / 4) * flux_val * (tau + 2. / 3);
    return bval;
}
