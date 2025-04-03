
typedef enum {EULER = 0, HEUN, RKFOUR, UNKNOWN} ode_t;

// ODE solver macros; __VA_ARGS__ are the parameters to the ODE fun

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
				for (size_t i=0; i<N; ++i) v[i] = u[i]+h2*udot1[i]; \
				odefun(udot2,v,N,__VA_ARGS__); \
				for (size_t i=0; i<N; ++i) v[i] = u[i]+h2*udot2[i]; \
				odefun(udot3,v,N,__VA_ARGS__); \
				for (size_t i=0; i<N; ++i) v[i] = u[i]+h*udot3[i]; \
				odefun(udot4,v,N,__VA_ARGS__); \
				double* const u1 = u+N; \
				for (size_t i=0; i<N; ++i) u1[i]  += u[i] + h6*(udot1[i]+2.0*udot2[i]+2.0*udot3[i]+udot4[i]); \
			}} \
			break; \
		default: \
			break; \
	} \
}
