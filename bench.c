#ifdef BENCH

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "stokes.h"
#include "advection_diffusion.h"


int 
main() 
{

    double current_time;

    double dx = (XMAX - XMIN) / (NX - 1);
    double dy = (YMAX - YMIN) / (NY - 1);

    //double x[NX];
    //double y[NY];
    double u[NY*NX];    // vel in x
    double v[NY*NX];    // vel in y
    double p[NY*NX];    // pressure
    double rho[NY*NX];  // density
    double nu[NY*NX];   // viscosity
    double t[NY*NX];    // temperature
    double k[NY*NX];    // conductivity

    double cp = 60.;
    double H = 0.;
    double dt = 9e-5;
    double ref_rho = 100.;

    int i, j;
    int timestep;
    
    // Initial conditions
    for ( j = 0; j < NY; j++ ){
       for ( i = 0; i < NX; i++ ){
          u[NX*j + i] = 0.; 
          v[NX*j + i] = 0.; 
          p[NX*j + i] = ref_rho * dy * GRAVITY * (NY - j);  // Approximate lithostatic pressure
          rho[NX*j + i] = ref_rho; 
          nu[NX*j + i] = 1.;
          // Make the temp field unstable 
          if ( j < 4 && i < (int)(NX/2) )
              t[NX*j + i] = 1000.;
          else if ( j > NY-5 && i > (int)(NX/2) )
              t[NX*j + i] = 0.;
          else
              t[NX*j + i] = 500.;
          k[NX*j + i] = 100.;
       }
    }

    // Init the boundaries too
    apply_vel_boundary_conditions(u, v);
    apply_thermal_boundary_conditions(t);

    timestep = 0;
    current_time = 0;
    
    // Main loop
    while ( timestep < 1000 ) {
        // Solve thermal stuff first, so the flow equations have some meat to start with
        solve_advection_diffusion(t, u, v, dx, dy, rho, dt, cp, k, H);
        update_rho(rho, t);
        update_nu(nu, t);
        update_k(k, t);
        solve_flow(u, v, dx, dy, p, rho, nu, dt);

        current_time += dt;
        timestep++;
    }

    return 0;
}
#endif
