#include "advection_diffusion.h"
#include <math.h>

void 
apply_thermal_boundary_conditions(double *t)
{
    int i, j;
   
    for ( j = 0; j < NY; j++ ){
        // Temp left wall, "freeslip" (temp doesn't escape)
        t[NX*j + 0] = t[NX*j + 1];
        
        // Temp right wall, "freeslip" (temp doesn't escape)
        t[NX*j + (NX-1)] = t[NX*j + (NX-2)];
    }

    for ( i = 0; i < NX; i++ ){
        // Temp bottom wall, 
        t[0 + i] = 1000.;

        // Temp bottom wall, 
        t[NX*(NY-1) + i] = 0.;
    }
}


void 
solve_advection_diffusion(double *t, 
                          double *u, 
                          double *v,
                          double dx, 
                          double dy,
                          double *rho, 
                          double dt,
                          double cp, 
                          double *k,
                          double H)
{
    int i,j; 
    double tn[NY*NX];
    double kx;
    double ky;

    // Pre-compute some basics
    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double twodx = 2. * dx;
    double twody = 2. * dy;

    for ( j = 0; j < NY; j++ ){
       for ( i = 0; i < NX; i++ ){
           tn[NX*j + i] = t[NX*j + i];
       }
    }
   
 
    for ( j = 1; j < NY-1; j++ ){
        for ( i = 1; i < NX-1; i++){
           kx = k[NX*j + i] * (tn[NX*j + (i+1)] - 2.*tn[NX*j + i] + tn[NX*j + (i-1)]) / dx2;
           ky = k[NX*j + i] * (tn[NX*(j+1) + i] - 2.*tn[NX*j + i] + tn[NX*(j-1) + i]) / dy2;

           t[NX*j + i] = tn[NX*j + i] + dt * ((H + kx + ky)/(rho[NX*j + i] * cp) \
                     - (u[NX*j + i] * ( (tn[NX*j + (i+1)] - tn[NX*j + (i-1)]) / twodx )) \
                     - (v[NX*j + i] * ( (tn[NX*(j+1) + i] - tn[NX*(j-1) + i]) / twody )) );
       }
    }

    apply_thermal_boundary_conditions(t);
}



void 
update_nu(double *nu, 
          double *t)
{
    /* Calculate temperature dependent viscosity. Based on 
     * Frank-Kamenetski formulation.
     * Viscosity gets lower with increase in temperature
     */

    int i, j;

    double ref_nu = 1.;
    double ref_temp = 500.;
    double theta = 1.5;

    for ( j = 0; j < NY; j++ ){
       for ( i = 0; i < NX; i++ ){
           nu[NX*j + i] = ref_nu * exp(-theta * ((t[NX*j + i] - ref_temp)/ref_temp));
       }
    }
}


void 
update_rho(double *rho, 
           double *t)
{
    /* Calculate density change due to thermal expansion.
     * Using boissinesq approximation (rho change does not
     * effect velocity field). Change is linear.
     */

    int i, j;

    double ref_rho = 100.;
    double ref_temp = 500.;
    double thermal_expansivity = 0.001;

    for ( j = 0; j < NY; j++ ){
       for ( i = 0; i < NX; i++ ){
           rho[NX*j + i] = ref_rho * (1. - (thermal_expansivity * (t[NX*j + i] - ref_temp)));
       }
    }
}


void 
update_k(double *k, 
         double *t)
{
    /* Calculate temperature dependent conductivity. 
     * Using linear relationship where higher temp gives
     * lower conductivity. This means slabs (cold drips)
     * will suck up heat faster. It's also relatively natural.
     */

    int i, j;

    double ref_k = 100.;
    double ref_temp = 500.;
    double thermal_factor = 0.001;

    for ( j = 0; j < NY; j++ ){
       for ( i = 0; i < NX; i++ ){
           k[NX*j + i] = ref_k * (1. - (thermal_factor * (t[NX*j + i] - ref_temp)));
       }
    }
}

