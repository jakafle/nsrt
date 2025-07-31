// program to calculate the weights in gaussian quadrature
#include <stdio.h>
#include <math.h>
#define EPS 3.0e-11
#define PI 3.141592654

void gauss_legendre(float x1, float x2, float x[], float w[], int n);

int main(void){
	float x1 = 1.0;
	float x2 = 6.0;
	float x[4];
	float w[4];
	int n = 4;

	gauss_legendre(x1,x2,x,w,n);
	printf("%f", w[1]);

}

void gauss_legendre(float x1, float x2, float x[], float w[], int n){

	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;

	// roots in symmetric interval
	m = (n+1)/2;
	xm = 0.5*(x2+x1);
	xl = 0.5*(x2-x1);

	// looping over the roots that we require
	for (i=1;i<=m;i++){

		z = cos(PI*(i-0.25)/(n+0.5));
		// Use Newton-Raphson method to find the root
		do {

			p1 = 1.0;
			p2 = 0.0;
			for (j=1;j<=n;j++){

				p3 = p2;
				p2 = p1;
				// the legendre polynomial
				p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}

			pp = n*(z*p1-p2)/(z*z - 1.0);
			z1 = z;
			z = z1-p1/pp;
		}
		while (fabs(z-z1) > EPS);
		
		x[i] = xm-xl*z;
		x[n+1-i] = xm+xl*z;
		w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
		w[n+1-i] = w[i];
	
	}

}
