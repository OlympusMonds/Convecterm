#include "stokes.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void 
solve_pressure_poisson(double (*p)[NX], 
                       double dx, 
                       double dy, 
                       double dt, 
                       double (*u)[NX], 
                       double (*v)[NX],
                       double (*rho)[NX])
{
    int i, j;
    int stepcount;
    
    double diff, pdif, pnt;
    double b[NY][NX];
    double pn[NY][NX];

    // Pre-compute
    double inv_dt = 1. / dt;
    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double twodx = 2. * dx;
    double twody = 2. * dy;

    // Pre-solve b term.
    for ( j = 1; j < (NY-1); j++ ){
        for ( i = 1; i < (NX-1); i++){
            b[j][i] = (( rho[j][i] * dx2 * dy2 ) / ( 2 * (dx2 + dy2))) * 
                      ( 
                          inv_dt * 
                          (
                           (( u[j][i+1] - u[j][i-1] ) / twodx ) +
                           (( v[j+1][i] - v[j-1][i] ) / twody )
                          ) -
                          pow( (( u[j][i+1] - u[j][i-1] ) / twodx), 2) -
                          2 *
                          ((( u[j+1][i] - u[j-1][i] ) / twody ) *
                           (( v[j][i+1] - v[j][i-1] ) / twodx )) -
                          pow( (( v[j+1][i] - v[j-1][i] ) / twody ), 2)
                      );
        }
    } 

    diff = 1.0;
    stepcount = 0;
    while ( 1 ) {
        if ( stepcount >= 1e6 ){
            printf(COLOR_RESET);
            printf("Unable to solve poisson: sc = %d, diff = %e\n", stepcount, diff);
            exit( 1 );
        }
        if ( diff < 1e-5 && stepcount > 4 )  // Once converged, break out of loop.
            break;

        for ( j = 0; j < NY; j++ ){
            for ( i = 0; i < NX; i++){
                pn[j][i] = p[j][i];
            }
        }

        for ( j = 1; j < NY-1; j++ ){
            for ( i = 1; i < NX-1; i++){
                p[j][i] = (( pn[j][i+1] + pn[j][i-1] ) * dy2 + ( pn[j+1][i] + pn[j-1][i] ) * dx2 ) /
                          ( 2 * (dx2 + dy2) ) - b[j][i];
            }
        }

        // Boundary conditions
        for ( j = 0; j < NY; j++ ){
            p[j][0] = p[j][1];
            p[j][NX-1] = p[j][NX-2];
        }
        for ( i = 0; i < NX; i++ ){
            p[0][i] = p[1][i];
            p[NY-1][i] = p[NY-2][i];
        }

        // Check if in steady state
        pdif = 0.;
        pnt = 0.;

        for ( j = 0; j < NY; j++ ){
           for ( i = 0; i < NX; i++ ){
               pdif += fabs(fabs(p[j][i]) - fabs(pn[j][i]));
               pnt += fabs(pn[j][i]);
           }
        }
        if ( pnt != 0. )
            diff = fabs(pdif/pnt);
        else
            diff = 1.0;

        stepcount += 1;
    }
}


void 
solve_stokes_momentum(double (*u)[NX], 
                      double (*v)[NX],
                      double (*un)[NX], 
                      double (*vn)[NX],
                      double (*p)[NX], 
                      double (*rho)[NX],
                      double (*nu)[NX],
                      double dt, 
                      double dx, 
                      double dy)
{
    int i, j;

    // Pre-compute
    double dtodx2 = dt / (dx * dx);
    double dtody2 = dt / (dy * dy);
   
    for ( j = 1; j < (NY-1); j++ ){
        for ( i = 1; i < (NX-1); i++){
            u[j][i] = un[j][i] - ( dt / (rho[j][i] * 2. * dx) ) * (p[j][i+1] - p[j][i-1]) +
                      nu[j][i] * (
                                  (dtodx2 * (un[j][i+1] - 2*un[j][i] + un[j][i-1])) +
                                  (dtody2 * (un[j+1][i] - 2*un[j][i] + un[j-1][i]))
                                 );
            
            v[j][i] = (vn[j][i] - ( dt / (rho[j][i] * 2. * dy) ) * (p[j+1][i] - p[j-1][i]) +
                       nu[j][i] * (
                                   (dtodx2 * (vn[j][i+1] - 2*vn[j][i] + vn[j][i-1])) +
                                   (dtody2 * (vn[j+1][i] - 2*vn[j][i] + vn[j-1][i]))
                                  )) - (GRAVITY * rho[j][i]);
        }
    }
}


void 
apply_vel_boundary_conditions(double (*u)[NX], 
                              double (*v)[NX])
{
    int i, j;
    
    for ( j = 0; j < NY; j++ ){
        // Vel Left wall, freeslip
        u[j][0] = 0.;
        v[j][0] = v[j][1];

        // Right wall, freeslip
        u[j][NX-1] = 0.;
        v[j][NX-1] = v[j][NX-2];
    }

    for ( i = 0; i < NX; i++ ){
        // Bottom wall, freeslip
        u[0][i] = u[1][i];
        v[0][i] = 0.;

        // Top wall, freeslip
        u[NY-1][i] = u[NY-2][i];
        v[NY-1][i] = 0.;
    }
}


void 
solve_flow(double (*u)[NX], 
           double (*v)[NX], 
           double dx, 
           double dy,
           double (*p)[NX], 
           double (*rho)[NX],
           double (*nu)[NX],
           double dt)
{

    int i, j;
    int stepcount = 0;

    double un[NY][NX];
    double vn[NY][NX];
    double diff = 1000.;
    double udif, vdif, unt, vnt;

    while ( 1 ) {
        if ( stepcount >= 50000 ){
            printf(COLOR_RESET);
            printf("Unable to solve: sc = %d, diff = %e\n", stepcount, diff);
            exit( 1 );
        }
        if ( diff < 1e-5 && stepcount >= 2 )
            break;

        for ( j = 0; j < NY; j++ ){
           for ( i = 0; i < NX; i++ ){
               un[j][i] = u[j][i];
               vn[j][i] = v[j][i];
           }
        }

        solve_pressure_poisson(p, dx, dy, dt, u, v, rho);
        solve_stokes_momentum(u, v, un, vn, p, rho, nu, dt, dx, dy);

        apply_vel_boundary_conditions(u, v);

        // Check if in steady state
        udif = 0.;
        vdif = 0.;
        unt = 0.;
        vnt = 0.;

        for ( j = 0; j < NY; j++ ){
           for ( i = 0; i < NX; i++ ){
               udif += fabs(fabs(u[j][i]) - fabs(un[j][i]));
               unt += fabs(un[j][i]);
               
               vdif += fabs(fabs(v[j][i]) - fabs(vn[j][i]));
               vnt += fabs(vn[j][i]);
           }
        }
        if ( unt != 0. && vnt != 0. )
            diff = fabs( fabs(udif/unt) - fabs(vdif/vnt) ) / 2.;
        else
            diff = 1.0;

        stepcount++;
    } 
}
