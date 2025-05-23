// odesolve test and demonstration program.
//
// Use Makefile in this directory to build.
//
// Results are plotted to file; if Gnuplot is available on your system, results are plotted;

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ode.h"
#include "mt64.h"

// The "Lorenz 96" chaotic system (https://en.wikipedia.org/wiki/Lorenz_96_model)

static inline void lorenz96(double* const xdot, const double* const x, const size_t N, const double F)
{
	xdot[0] = (x[1]-x[N-2])*x[N-1]-x[0]+F;
	xdot[1] = (x[2]-x[N-1])*x[0]-x[1]+F;
	for (size_t i=2; i<N-1; ++i) xdot[i] = (x[i+1]-x[i-2])*x[i-1]-x[i]+F;
	xdot[N-1] = (x[0]-x[N-3])*x[N-2]-x[N-1]+F;
}

// Main function

int lorenz96test(int argc, char* argv[])
{
	// Default command-line parameters

	const double      F   = argc > 1 ?         atof(argv[1])   : 8.0;    // Lorenz 96 F parameter
	const size_t      N   = argc > 2 ? (size_t)atol(argv[2])   : 5;      // system dimension (number of variables)
	const double      dt  = argc > 3 ?         atof(argv[3])   : 0.01;   // integration time step
	const size_t      n   = argc > 4 ? (size_t)atol(argv[4])   : 10000;  // number of integration time steps
	const char* const ode = argc > 5 ?              argv[5]    : "Heun"; // "Euler", "Heun", or "RK4"
	const char* const of  = argc > 6 ?              argv[6]    : "/tmp/lorenz96.asc";
#ifdef HAVE_GNUPLOT
	const char* const gf  = argc > 7 ?              argv[7]    : "/tmp/lorenz96.gp";
#endif //HAVE_GNUPLOT

	// Display command-line  parameters

	printf("\n*** ODESOLVE test (Lorenz 96 system) ***\n\n");
	printf("system dimension            =  %zu\n",  N);
	printf("Lorenz 96 F parameter       =  %g\n",   F);
	printf("integration step size       =  %g\n",   dt);
	printf("number of integration steps =  %zu\n",  n);
	printf("ODE solver                  =  %s\n\n", ode);

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

// Ornstein-Uhlenbeck process (stochastic differential equation: https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process)

static inline double ouproc(const double x, const double a)
{
	return -a*x;
}

int outest(int argc, char* argv[])
{
	// Default command-line parameters

	const double      a    = argc > 1 ?           atof(argv[1])   : 0.1;    // OU decay parameter
	const double      sig  = argc > 2 ?           atof(argv[2])   : 1.0;    // OU Wiener noise intensity
	const double      dt   = argc > 3 ?           atof(argv[3])   : 0.01;   // integration time step
	const size_t      n    = argc > 4 ?   (size_t)atol(argv[4])   : 10000;  // number of integration time steps
	const mtuint_t    seed = argc > 5 ? (mtuint_t)atol(argv[5])   : 0;      // PRNG seed (0 for random random seed :-)
	const char* const ode  = argc > 6 ?                argv[6]    : "Heun"; // "Euler", "Heun", or "RK4"
	const char* const of   = argc > 7 ?                argv[7]    : "/tmp/ou.asc";
#ifdef HAVE_GNUPLOT
	const char* const gf   = argc > 8 ?                argv[8]    : "/tmp/ou.gp";
#endif //HAVE_GNUPLOT

	// Display command-line  parameters

	printf("\n*** ODESOLVE test (Ornstein-Uhlenbeck process) ***\n\n");
	printf("OU decay parameter          = %g\n",    a   );
	printf("OU noise intensity          = %g\n",    sig );
	printf("integration step size       = %g\n",    dt  );
	printf("number of integration steps = %zu\n",   n   );
	printf("random seed                 = %zu%seed\n",seed,seed?"":" (random random seed :-)");
	printf("ODE solver                  = %s\n\n",  ode );

	// Check command-line parameters

	const ode_t solver = str2ode(ode);
	if (solver == UNKNOWN) {
		fprintf(stderr,"ERROR: Unknown ODE solver\n");
		return EXIT_FAILURE;
	}

	// Mersenne Twister pseudo-random number generator

	mt_t rng;           // the PRNG
	mt_seed(&rng,seed); // initialise the PRNG

	// Allocate memory for OU variable

	double* const x = calloc(n,sizeof(double));

	// Prefill variable with Wiener noise (note that *variance* -- not std. dev.! -- scales linearly with time increment)

	const double ssig = sig*sqrt(dt); // scaled noise std. dev.
	for (size_t i=1; i<n; ++i) x[i] = ssig*mt_randn(&rng);

	// Solve the ODE

	ODE1(solver,ouproc,x,n,dt,a);

	// Write results to file

	FILE* const offs = fopen(of,"w");
	if (offs == NULL) {
		perror("ERROR: Failed to open output file");
		return EXIT_FAILURE;
	}
	for (size_t i=0; i<n; ++i) {
		fprintf(offs,"%16.8f %16.8f\n",((double)(i+1))*dt, x[i]);
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
	fprintf(gpfs,"set title \"Ornstein-Uhlenbeck process (%s solver)\"\n",ode);
	fprintf(gpfs,"set xlabel \"t (time)\"\n");
	fprintf(gpfs,"set ylabel \"x\"\n");
	fprintf(gpfs,"plot \"%s\" u 1:2 w l not\n",of);
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

// Main function

static const int ntests = 2;

int main(int argc, char* argv[])
{
	if (argc < 2) {
		fprintf(stderr,"\n*** Must specify test (1 - %d) ***\n\n",ntests);
		return EXIT_FAILURE;
	}
	const int test = atoi(argv[1]);
	if (test < 1 || test > ntests) {
		fprintf(stderr,"Test number must be 1 - %d\n",ntests);
		return EXIT_FAILURE;
	}

	switch (test) {
		case 1 : return lorenz96test (argc-1,argv+1);
		case 2 : return outest       (argc-1,argv+1);
	}
	return EXIT_FAILURE; // shouldn't get here!
}
