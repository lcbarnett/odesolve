#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ode.h"

#define UNUSED __attribute__ ((unused))

// The "Lorenz 96" chaotic system

inline void lorenz96(double* const xdot, const double* const x, const size_t N, const double F)
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

	const double      F     = argc > 1 ?         atof(argv[1])   : 8.0;
	const size_t      N     = argc > 2 ? (size_t)atol(argv[2])   : 5;
	const double      dt    = argc > 3 ?         atof(argv[3])   : 0.01;
	const size_t      n     = argc > 4 ? (size_t)atol(argv[4])   : 10000;
	const char* const sstr  = argc > 5 ?              argv[5]    : "Heun";
	const char* const ofile = argc > 6 ?              argv[6]    : "/tmp/test.asc";
#ifdef HAVE_GNUPLOT
	const char* const gfile = argc > 6 ?              argv[6]    : "/tmp/test.gp";
#endif //HAVE_GNUPLOT

	if (N < 3)  {
		fprintf(stderr,"ERROR: Need at least three variables\n");
		return EXIT_FAILURE;
	}

	// Check the solver

	const ode_t solver =
		strcasecmp(sstr,"Euler") == 0 ? EULER  :
		strcasecmp(sstr,"Heun" ) == 0 ? HEUN   :
		strcasecmp(sstr,"RK4"  ) == 0 ? RKFOUR : UNKNOWN;

	if (solver == UNKNOWN) {
		fprintf(stderr,"ERROR: Unknown ODE solver\n");
		return EXIT_FAILURE;
	}

	// Display parameters

	printf("\n*** ODESOLVE test (Lorenz 96 system) ***\n\n");
	printf("Lorenz 96 dimension         =  %zu\n",   N);
	printf("Lorenz 96 F parameter       =  %g\n",    F);
	printf("integration step size       =  %g\n",    dt);
	printf("number of integration steps =  %zu\n",   n);
	printf("ODE solver                  =  %s\n\n",  sstr);

	// Allocate memory for variables

	double* const x = calloc(N*n,sizeof(double));

	// Set some initial values (here we need at least one variable not to be zero)

	x[0] = 1.0;

	// Solve the ODE

	ODE(solver,lorenz96,x,N,n,dt,F);

	// Write results to file

	FILE* const fp = fopen(ofile,"w");
	if (fp == NULL) {
		perror("ERROR: Failed to open output file");
		return EXIT_FAILURE;
	}
	for (size_t k=0; k<n; ++k) {
		const double* const xk = x + N*k;
		for (size_t i=0; i<N; ++i) fprintf(fp,"%17.8f",xk[i]);
		fputc('\n',fp);
	}
	if (fclose(fp) != 0) {
		perror("ERROR: Failed to close output file");
		return EXIT_FAILURE;
	}

	free(x); // finished with it

	// if Gnuplot available, plot trajectory of first three variables in 3D

#ifdef HAVE_GNUPLOT
	FILE* const gp = fopen(gfile,"w");
	if (gp == NULL) {
		perror("ERROR: failed to open Gnuplot command file\n");
		return EXIT_FAILURE;
	}
	fprintf(gp,"unset key\n");
	fprintf(gp,"set grid\n");
	fprintf(gp,"set title \"Lorenz 96 system (%s solver)\"\n",sstr);
	fprintf(gp,"set xlabel \"x\"\n");
	fprintf(gp,"set ylabel \"y\"\n");
	fprintf(gp,"set zlabel \"z\"\n");
	fprintf(gp,"splot \"%s\" u 1:2:3 w l not\n",ofile);
	if (fclose(gp) != 0) {
		perror("Failed to close Gnuplot command file");
		return EXIT_FAILURE;
	}
	const size_t strlen = 100;
	char gpcmd[strlen+1];
	snprintf(gpcmd,strlen,"gnuplot -p %s",gfile);
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
