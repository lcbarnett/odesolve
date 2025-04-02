#ifndef ODE_H
#define ODE_H

#include <math.h> // for maths functions

#define UNUSED __attribute__ ((unused))

// The "Lorenz 96" chaotic system

static inline void lorenz96(double* const xdot, const double* const x, const size_t N, const double F)
{
	xdot[0] = (x[1]-x[N-2])*x[N-1]-x[0]+F;
	xdot[1] = (x[2]-x[N-1])*x[0]-x[1]+F;
	for (size_t i=2; i<N-1; ++i) xdot[i] = (x[i+1]-x[i-2])*x[i-1]-x[i]+F;
	xdot[N-1] = (x[0]-x[N-3])*x[N-2]-x[N-1]+F;
}

// Main function

int main(int argc, char* argv[])
{
	// Default command-line parameters

	const double F      = argc > 1 ?         atof(argv[1])   : 8.0;
	const size_t N      = argc > 2 ? (size_t)atol(argv[2])   : 10000;
	const double dt     = argc > 3 ?         atof(argv[3])   : 0.001;
	const size_t n      = argc > 4 ? (size_t)atol(argv[4])   : 10000;
	const ode_t  solver = argc > 5 ?  (ode_t)atoi(argv[5])   : HEUN;

	if (solver < 0 || solver > RKFOUR) {
		fprintf(stderr,"Unknown ODE solver\n");
		return EXIT_FAILURE;
	}

	printf("\n*** ODESOLVE test (Lorenz 96 system) ***\n\n");
	printf("Lorenz 96 dimension         =  %zu\n",   N);
	printf("Lorenz 96 F parameter       =  %g\n",    F);
	printf("integration step size       =  %g\n",    h);
	printf("number of integration steps =  %zu\n",   n);
	printf("ODE solver                  =  %s\n",    ode2str(s));

	double* const x = calloc(N*n,sizeof(double));

	ODE(solver,lorenz96,x,N,n,h,F);

	return(EXIT_SUCCESS);
}

#endif // ODE_H
