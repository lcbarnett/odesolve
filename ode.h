// Generic ODE solver macro, with Euler, Heun and Runge-Kutta 4 ("RK4") integration.
//
// See test/test.c and test/Makefile for example usage, and for building programs using ode.h

#include <string.h>

typedef enum {EULER = 0, HEUN, RKFOUR, UNKNOWN} ode_t;

static inline ode_t str2ode(const char* const str)
{
	return strcasecmp(str,"Euler") == 0 ? EULER : strcasecmp(str,"Heun" ) == 0 ? HEUN : strcasecmp(str,"RK4"  ) == 0 ? RKFOUR : UNKNOWN;
}

// the ODE Macro parameters:
//
// Name     Description              Type
// ----------------------------------------------------------------------------------------------------------------
// ode      ODE type                 ode_t
// odefun   ODE function pointer     void (*odefun)(double* const xdot, const double* const x, const size_t N, ...)
// x        ODE variables            double* const
// N        Sytem dimension          const size_t
// n        Number of time steps     const size_t
// h        Integration step size    const double
// ...      'odefun' parameters      as specified in the 'odefun' prototype
// ----------------------------------------------------------------------------------------------------------------
//
// Memory for the ODE variables should be allocated with
//
//   	double* const x = calloc(N*n,sizeof(double));
//
// then initialised appropriately, and deallocated after use with free(x).

#define ODE(ode,odefun,x,N,n,h,...) \
{ \
	switch (ode) { \
		case EULER: { \
			printf("EULER : "#odefun"\n"); \
			double udot[N]; \
			for (double* u=x; u<x+N*(n-1); u+=N) { \
				odefun(udot,u,N,__VA_ARGS__); \
				double* const u1 = u+N; \
				for (size_t i=0; i<N; ++i) u1[i] += u[i] + h*udot[i]; \
			}} \
			break; \
		case HEUN: { \
			printf("HEUN : "#odefun"\n"); \
			const double h2 = h/2.0; \
			double udot1[N], udot2[N]; \
			double v[N]; \
			for (double* u=x; u<x+N*(n-1); u+=N) { \
				odefun(udot1,u,N,__VA_ARGS__); \
				for (size_t i=0; i<N; ++i) v[i] = u[i] + h*udot1[i]; \
				odefun(udot2,v,N,__VA_ARGS__); \
				double* const u1 = u+N; \
				for (size_t i=0; i<N; ++i) u1[i] += u[i] + h2*(udot1[i]+udot2[i]); \
			}} \
			break; \
		case RKFOUR: { \
			printf("RK4 : "#odefun"\n"); \
			const double h2 = h/2.0; \
			const double h6 = h/6.0; \
			double udot1[N],udot2[N],udot3[N],udot4[N]; \
			double v[N]; \
			for (double* u=x; u<x+N*(n-1); u+=N) { \
				odefun(udot1,u,N,__VA_ARGS__); \
				for (size_t i=0; i<N; ++i) v[i] = u[i] + h2*udot1[i]; \
				odefun(udot2,v,N,__VA_ARGS__); \
				for (size_t i=0; i<N; ++i) v[i] = u[i] + h2*udot2[i]; \
				odefun(udot3,v,N,__VA_ARGS__); \
				for (size_t i=0; i<N; ++i) v[i] = u[i] + h*udot3[i]; \
				odefun(udot4,v,N,__VA_ARGS__); \
				double* const u1 = u+N; \
				for (size_t i=0; i<N; ++i) u1[i]  += u[i] + h6*(udot1[i]+2.0*udot2[i]+2.0*udot3[i]+udot4[i]); \
			}} \
			break; \
		default: \
			break; \
	} \
}

// More efficient for 1-dimesnional ODEs (N = 1)

// the ODE1 Macro parameters; same, except no N parameter, and
//
// Name     Description              Type
// ----------------------------------------------------------------------------------------------------------------
// odefun   ODE function pointer     void (*odefun)(double* const xdot, const double* const x, ...)
// ----------------------------------------------------------------------------------------------------------------

#define ODE1(ode,odefun,x,n,h,...) \
{ \
	switch (ode) { \
		case EULER: { \
			printf("EULER : "#odefun"\n"); \
			double udot; \
			for (double* u=x; u<x+n-1; ++u) { \
				odefun(udot,u,__VA_ARGS__); \
				*(u+1) += *u + h*udot; \
			}} \
			break; \
		case HEUN: { \
			printf("HEUN : "#odefun"\n"); \
			const double h2 = h/2.0; \
			double udot1, udot2; \
			double v; \
			for (double* u=x; u<x+n-1; ++u) { \
				odefun(udot1,u,__VA_ARGS__); \
				v = *u + h*udot1; \
				odefun(udot2,v,__VA_ARGS__); \
				*(u+1) += *u + h2*(udot1+udot2); \
			}} \
			break; \
		case RKFOUR: { \
			printf("RK4 : "#odefun"\n"); \
			const double h2 = h/2.0; \
			const double h6 = h/6.0; \
			double udot1,udot2,udot3,udot4; \
			double v; \
			for (double* u=x; u<x+n-1; ++u) { \
				odefun(udot1,u,__VA_ARGS__); \
				v = *u + h2*udot1; \
				odefun(udot2,v,__VA_ARGS__); \
				v = *u + h2*udot2; \
				odefun(udot3,v,__VA_ARGS__); \
				v = *u + h*udot3; \
				odefun(udot4,v,__VA_ARGS__); \
				*(u+1) += u + h6*(udot1+2.0*udot2+2.0*udot3+udot4); \
			}} \
			break; \
		default: \
			break; \
	} \
}
