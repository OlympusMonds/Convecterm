#include "advection_diffusion.h"
#include <math.h>

void 
apply_thermal_boundary_conditions(double (*t)[NX])
{
    int i, j;
    
    for ( j = 0; j < NY; j++ ){
        // Temp left wall, "freeslip" (temp doesn't escape)
        t[j][0] = t[j][1];
        
        // Temp right wall, "freeslip" (temp doesn't escape)
        t[j][NX-1] = t[j][NX-2];
    }

    for ( i = 0; i < NX; i++ ){
        // Temp bottom wall, 
        t[0][i] = 1000.;

        // Temp bottom wall, 
        t[NY-1][i] = 0.;
    }
}


void 
solve_advection_diffusion(double (*t)[NX], 
                          double (*u)[NX], 
                          double (*v)[NX],
                          double dx, 
                          double dy,
                          double (*rho)[NX], 
                          double dt,
                          double cp, 
                          double (*k)[NX],
                          double H)
{
    int i,j; 
    double tn[NY][NX];
    double kx;
    double ky;

    // Pre-compute some basics
    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double twodx = 2. * dx;
    double twody = 2. * dy;

    for ( j = 0; j < NY; j++ ){
       for ( i = 0; i < NX; i++ ){
           tn[j][i] = t[j][i];
       }
    }
   
 
    for ( j = 1; j < NY-1; j++ ){
        for ( i = 1; i < NX-1; i++){
           kx = k[j][i] * (tn[j][i+1] - 2.*tn[j][i] + tn[j][i-1]) / dx2;
           ky = k[j][i] * (tn[j+1][i] - 2.*tn[j][i] + tn[j-1][i]) / dy2;

           t[j][i] = tn[j][i] + dt * ((H + kx + ky)/(rho[j][i] * cp) \
                     - (u[j][i] * ( (tn[j][i+1] - tn[j][i-1]) / twodx )) \
                     - (v[j][i] * ( (tn[j+1][i] - tn[j-1][i]) / twody )) );
       }
    }

    apply_thermal_boundary_conditions(t);
}



void 
update_nu(double (*nu)[NX], 
          double (*t)[NX])
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
           nu[j][i] = ref_nu * exp(-theta * ((t[j][i] - ref_temp)/ref_temp));
       }
    }
}


void 
update_rho(double (*rho)[NX], 
           double (*t)[NX])
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
           rho[j][i] = ref_rho * (1. - (thermal_expansivity * (t[j][i] - ref_temp)));
       }
    }
}


void 
update_k(double (*k)[NX], 
         double (*t)[NX])
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
           k[j][i] = ref_k * (1. - (thermal_factor * (t[j][i] - ref_temp)));
       }
    }
}

