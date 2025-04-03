// odesolve test and demonstration program.
//
// Use Makefile in this directory to build.
//
// Results are plotted to file; if Gnuplot is available on your system, results are plotted;

// The "Lorenz 96" chaotic system (https://en.wikipedia.org/wiki/Lorenz_96_model)

#include <stdlib.h>
#include <stdio.h>

#include "ode.h"

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

	const double      F   = argc > 1 ?         atof(argv[1])   : 8.0;    // Lorenz 96 F parameter
	const size_t      N   = argc > 2 ? (size_t)atol(argv[2])   : 5;      // system dimension (number of variables)
	const double      dt  = argc > 3 ?         atof(argv[3])   : 0.01;   // integration time step
	const size_t      n   = argc > 4 ? (size_t)atol(argv[4])   : 10000;  // number of integration time steps
	const char* const ode = argc > 5 ?              argv[5]    : "Heun"; // "Euler", "Heun", or "RK4"
	const char* const of  = argc > 6 ?              argv[6]    : "/tmp/test.asc";
#ifdef HAVE_GNUPLOT
	const char* const gf  = argc > 7 ?              argv[7]    : "/tmp/test.gpfs";
#endif //HAVE_GNUPLOT

	// Display command-line  parameters

	printf("\n*** ODESOLVE test (Lorenz 96 system) ***\n\n");
	printf("system dimension            =  %zu\n",   N);
	printf("Lorenz 96 F parameter       =  %g\n",    F);
	printf("integration step size       =  %g\n",    dt);
	printf("number of integration steps =  %zu\n",   n);
	printf("ODE solver                  =  %s\n\n",  ode);

	// Check command-line parameters

	if (N < 4)  {
		fprintf(stderr,"ERROR: Lorenz 96 needs at least four variables\n");
		return EXIT_FAILURE;
	}

	const ode_t solver = str2ode(ode);
	if (solver == UNKNOWN) {
		fprintf(stderr,"ERROR: Unknown ODE solver\n");
		return EXIT_FAILURE;
	}

	// Allocate memory for variables

	double* const x = calloc(N*n,sizeof(double));

	// Set some initial values (here we need at least one variable not to be zero)

	x[0] = 1.0;

	// Solve the ODE

	ODE(solver,lorenz96,x,N,n,dt,F);

	// Write results to file

	FILE* const offs = fopen(of,"w");
	if (offs == NULL) {
		perror("ERROR: Failed to open output file");
		return EXIT_FAILURE;
	}
	for (size_t k=0; k<n; ++k) {
		const double* const xk = x + N*k;
		for (size_t i=0; i<N; ++i) fprintf(offs," %16.8f",xk[i]);
		fputc('\n',offs);
	}
	if (fclose(offs) != 0) {
		perror("ERROR: Failed to close output file");
		return EXIT_FAILURE;
	}

	free(x); // finished with it

	// if Gnuplot available, plot trajectory of first three variables in 3D

#ifdef HAVE_GNUPLOT
	FILE* const gpfs = fopen(gf,"w");
	if (gpfs == NULL) {
		perror("ERROR: failed to open Gnuplot command file\n");
		return EXIT_FAILURE;
	}
	fprintf(gpfs,"unset key\n");
	fprintf(gpfs,"set grid\n");
	fprintf(gpfs,"set title \"Lorenz 96 system (%s solver)\"\n",ode);
	fprintf(gpfs,"set xlabel \"x\"\n");
	fprintf(gpfs,"set ylabel \"y\"\n");
	fprintf(gpfs,"set zlabel \"z\"\n");
	fprintf(gpfs,"splot \"%s\" u 1:2:3 w l not\n",of);
	if (fclose(gpfs) != 0) {
		perror("Failed to close Gnuplot command file");
		return EXIT_FAILURE;
	}
	const size_t strlen = 100;
	char gpcmd[strlen+1];
	snprintf(gpcmd,strlen,"gnuplot -p %s",gf);
	printf("\nGnuplot command: %s\n\n",gpcmd);
	if (system(gpcmd) == -1) {
		perror("ERROR: Failed to run Gnuplot command");
		return EXIT_FAILURE;
	}
#else
	printf("\nNOTE: Gnuplot unavailable: can't plot\n\n");
#endif //HAVE_GNUPLOT

	return EXIT_SUCCESS;
}
