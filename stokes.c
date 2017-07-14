#include "stokes.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void 
solve_pressure_poisson(double *p, 
                       double dx, 
                       double dy, 
                       double dt, 
                       double *u, 
                       double *v,
                       double *rho)
{
    int i, j;
    int stepcount;
    
    double diff, pdif, pnt;
    double b[NY*NX];
    double pn[NY*NX];

    // Pre-compute
    double inv_dt = 1. / dt;
    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double twodx = 2. * dx;
    double twody = 2. * dy;

    // Pre-solve b term.
    #pragma omp parallel for
    for ( j = 1; j < (NY-1); j++ ){
        for ( i = 1; i < (NX-1); i++){
            b[NX*j + i] = (( rho[NX*j + i] * dx2 * dy2 ) / ( 2 * (dx2 + dy2))) * 
                      ( 
                          inv_dt * 
                          (
                           (( u[NX*j + (i+1)] - u[NX*j + (i-1)] ) / twodx ) +
                           (( v[NX*(j+1) + i] - v[NX*(j-1) + i] ) / twody )
                          ) -
                          pow( (( u[NX*j + (i+1)] - u[NX*j + (i-1)] ) / twodx), 2) -
                          2 *
                          ((( u[NX*(j+1) + i] - u[NX*(j-1) + i] ) / twody ) *
                           (( v[NX*j + (i+1)] - v[NX*j + (i-1)] ) / twodx )) -
                          pow( (( v[NX*(j+1) + i] - v[NX*(j-1) + i] ) / twody ), 2)
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

        #pragma omp parallel for
        for ( j = 0; j < NY; j++ ){
            for ( i = 0; i < NX; i++){
                pn[NX*j + i] = p[NX*j + i];
            }
        }

        #pragma omp parallel for
        for ( j = 1; j < NY-1; j++ ){
            for ( i = 1; i < NX-1; i++){
                p[NX*j + i] = (( pn[NX*j + (i+1)] + pn[NX*j + (i-1)] ) * dy2 + ( pn[NX*(j+1) + i] + pn[NX*(j-1) + i] ) * dx2 ) /
                          ( 2 * (dx2 + dy2) ) - b[NX*j + i];
            }
        }

        // Boundary conditions
        #pragma omp parallel for
        for ( j = 0; j < NY; j++ ){
            p[NX*j + 0] = p[NX*j + 1];
            p[NX*j + (NX-1)] = p[NX*j + (NX-2)];
        }
        #pragma omp parallel for
        for ( i = 0; i < NX; i++ ){
            p[0 + i] = p[NX*1 + i];
            p[NX*(NY-1) + i]= p[NX*(NY-2) + i];
        }

        // Check if in steady state
        pdif = 0.;
        pnt = 0.;

        #pragma omp parallel for
        for ( j = 0; j < NY; j++ ){
           for ( i = 0; i < NX; i++ ){
               pdif += fabs(fabs(p[NX*j + i]) - fabs(pn[NX*j + i]));
               pnt += fabs(pn[NX*j + i]);
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
solve_stokes_momentum(double *u, 
                      double *v,
                      double *un, 
                      double *vn,
                      double *p, 
                      double *rho,
                      double *nu,
                      double dt, 
                      double dx, 
                      double dy)
{
    int i, j;

    // Pre-compute
    double dtodx2 = dt / (dx * dx);
    double dtody2 = dt / (dy * dy);
   
    #pragma omp parallel for
    for ( j = 1; j < (NY-1); j++ ){
        for ( i = 1; i < (NX-1); i++){
            u[NX*j + i] = un[NX*j + i] - ( dt / (rho[NX*j + i] * 2. * dx) ) * (p[NX*j + (i+1)] - p[NX*j + (i-1)]) +
                      nu[NX*j + i] * (
                                  (dtodx2 * (un[NX*j + (i+1)] - 2*un[NX*j + i] + un[NX*j + (i-1)])) +
                                  (dtody2 * (un[NX*(j+1) + i] - 2*un[NX*j + i] + un[NX*(j-1) + i]))
                                 );
            
            v[NX*j + i] = (vn[NX*j + i] - ( dt / (rho[NX*j + i] * 2. * dy) ) * (p[NX*(j+1) + i] - p[NX*(j-1) + i]) +
                       nu[NX*j + i] * (
                                   (dtodx2 * (vn[NX*j + (i+1)] - 2*vn[NX*j + i] + vn[NX*j + (i-1)])) +
                                   (dtody2 * (vn[NX*(j+1) + i] - 2*vn[NX*j + i] + vn[NX*(j-1) + i]))
                                  )) - (GRAVITY * rho[NX*j + i]);
        }
    }
}


void 
apply_vel_boundary_conditions(double *u, 
                              double *v)
{
    int i, j;
    
    #pragma omp parallel for
    for ( j = 0; j < NY; j++ ){
        // Vel Left wall, freeslip
        u[NX*j + 0] = 0.;
        v[NX*j + 0] = v[NX*j + 1];

        // Right wall, freeslip
        u[NX*j + (NX-1)] = 0.;
        v[NX*j + (NX-1)] = v[NX*j + (NX-2)];
    }

    #pragma omp parallel for
    for ( i = 0; i < NX; i++ ){
        // Bottom wall, freeslip
        u[0 + i] = u[NX*1 + i];
        v[0 + i] = 0.;

        // Top wall, freeslip
        u[NX*(NY-1) + i] = u[NX*(NY-2) + i];
        v[NX*(NY-1) + i] = 0.;
    }
}


void 
solve_flow(double *u, 
           double *v, 
           double dx, 
           double dy,
           double *p, 
           double *rho,
           double *nu,
           double dt)
{

    int i, j;
    int stepcount = 0;

    double un[NY*NX];
    double vn[NY*NX];
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

        #pragma omp parallel for
        for ( j = 0; j < NY; j++ ){
           for ( i = 0; i < NX; i++ ){
               un[NX*j + i] = u[NX*j + i];
               vn[NX*j + i] = v[NX*j + i];
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

        #pragma omp parallel for
        for ( j = 0; j < NY; j++ ){
           for ( i = 0; i < NX; i++ ){
               udif += fabs(fabs(u[NX*j + i]) - fabs(un[NX*j + i]));
               unt += fabs(un[NX*j + i]);
               
               vdif += fabs(fabs(v[NX*j + i]) - fabs(vn[NX*j + i]));
               vnt += fabs(vn[NX*j + i]);
           }
        }
        if ( unt != 0. && vnt != 0. )
            diff = fabs( fabs(udif/unt) - fabs(vdif/vnt) ) / 2.;
        else
            diff = 1.0;

        stepcount++;
    } 
}
