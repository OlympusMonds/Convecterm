#include "stokes.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void 
solve_pressure_poisson(double *restrict p, 
                       double dx, 
                       double dy, 
                       double dt, 
                       double *restrict u, 
                       double *restrict v,
                       double *restrict rho)
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
    unsigned int n = 0;
    unsigned int s = 0;
    unsigned int m = 0;
    unsigned int e = 0;
    unsigned int w = 0;


    // Pre-solve b term.
    #pragma omp parallel for
    for ( j = 1; j < (NY-1); j++ ){
        #pragma omp simd safelen(3)
        for ( i = 1; i < (NX-1); i++){
            m = NX*j + i;
            n = NX*(j-1) + i;
            s = NX*(j+1) + i;
            e = NX*j + (i+1);
            w = NX*j + (i-1);

            b[m] = (( rho[m] * dx2 * dy2 ) / ( 2 * (dx2 + dy2))) * 
                      ( 
                          inv_dt * 
                          (
                           (( u[e] - u[w] ) / twodx ) +
                           (( v[s] - v[n] ) / twody )
                          ) -
                          pow( (( u[e] - u[w] ) / twodx), 2) -
                          2 *
                          ((( u[s] - u[n] ) / twody ) *
                           (( v[e] - v[w] ) / twodx )) -
                          pow( (( v[s] - v[n] ) / twody ), 2)
                      );
        }
    } 

    diff = 1.0;
    stepcount = 0;
    while ( 1 ) {
        if ( stepcount >= 1e7 ){
            printf(COLOR_RESET);
            printf("Unable to solve poisson: sc = %d, diff = %e\n", stepcount, diff);
            exit( 1 );
        }
        if ( diff < 1e-5 && stepcount > 4 )  // Once converged, break out of loop.
            break;

        #pragma omp parallel for
        for ( i = 0; i < NX*NY; i++){
            pn[i] = p[i];
        }

        #pragma omp parallel for
        for ( j = 1; j < NY-1; j++ ){
            #pragma omp simd safelen(3)
            for ( i = 1; i < NX-1; i++){
                m = NX*j + i;
                n = NX*(j-1) + i;
                s = NX*(j+1) + i;
                e = NX*j + (i+1);
                w = NX*j + (i-1);
                
                p[m] = (( pn[e] + pn[w] ) * dy2 + ( pn[s] + pn[n] ) * dx2 ) /
                          ( 2 * (dx2 + dy2) ) - b[m];
            }
        }

        // Boundary conditions
        #pragma omp parallel for
        for ( j = 0; j < NY; j++ ){
            p[NX*j + 0] = p[NX*j + 1];
            p[NX*j + (NX-1)] = p[NX*j + (NX-2)];
        }
        // Periodic boundary condition for vert walls
        /*
        #pragma omp parallel for
        for ( j = 0; j < NY; j++ ){
            p[NX*j + 0] = p[NX*j + (NX-1)];
        }
        */
        #pragma omp parallel for
        for ( i = 0; i < NX; i++ ){
            //p[0 + i] = p[NX*1 + i];
            //p[NX*(NY-1) + i]= p[NX*(NY-2) + i];
            p[0 + i] = 0.;
            p[NX*(NY-1) + i] = 0.;
        }

        // Check if in steady state
        pdif = 0.;
        pnt = 0.;

        #pragma omp parallel for
        for ( i = 0; i < NX*NY; i++ ){
            pdif += fabs(fabs(p[i]) - fabs(pn[i]));
            pnt += fabs(pn[i]);
        }

        if ( pnt > 0. )
            diff = fabs(pdif/pnt);
        else
            diff = 1.0;

        if ( stepcount % 10000 == 0) {
            printf("sc: %d, pdif: %e, pnt: %e, Diff: %e\n", stepcount, pdif, pnt, diff);
        }
        stepcount += 1;
    }
}


void 
solve_stokes_momentum(double *restrict u, 
                      double *restrict v,
                      double *restrict un, 
                      double *restrict vn,
                      double *restrict p, 
                      double *restrict rho,
                      double *restrict nu,
                      double dt, 
                      double dx, 
                      double dy)
{
    int i, j;

    unsigned int n = 0;
    unsigned int s = 0;
    unsigned int m = 0;
    unsigned int e = 0;
    unsigned int w = 0;


    // Pre-compute
    double dtodx2 = dt / (dx * dx);
    double dtody2 = dt / (dy * dy);
   
    #pragma omp parallel for
    for ( j = 1; j < (NY-1); j++ ){
        #pragma omp simd safelen(3)
        for ( i = 1; i < (NX-1); i++){
            m = NX*j + i;
            n = NX*(j-1) + i;
            s = NX*(j+1) + i;
            e = NX*j + (i+1);
            w = NX*j + (i-1);

            u[m] = un[m] - ( dt / (rho[m] * 2. * dx) ) * (p[e] - p[w]) +
                      nu[m] * (
                                  (dtodx2 * (un[e] - 2*un[m] + un[w])) +
                                  (dtody2 * (un[s] - 2*un[m] + un[n]))
                                 );
            
            v[m] = (vn[m] - ( dt / (rho[m] * 2. * dy) ) * (p[s] - p[n]) +
                       nu[m] * (
                                   (dtodx2 * (vn[e] - 2*vn[m] + vn[w])) +
                                   (dtody2 * (vn[s] - 2*vn[m] + vn[n]))
                                  )) - (GRAVITY * rho[m]);
        }
    }
}


void 
apply_vel_boundary_conditions(double *restrict u, 
                              double *restrict v)
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
    /*  
    #pragma omp parallel for
    for ( j = 0; j < NY; j++ ){
        // Periodic BCs
        u[NX*j + 0] = u[NX*j + (NX-1)];
        v[NX*j + 0] = v[NX*j + (NX-1)];
    }
    */

    #pragma omp parallel for
    for ( i = 0; i < NX; i++ ){
        // Bottom wall, freeslip
        u[0 + i] = 0.; //u[NX*1 + i];
        v[0 + i] = 0.;

        // Top wall, freeslip
        u[NX*(NY-1) + i] = 0.; //u[NX*(NY-2) + i];
        v[NX*(NY-1) + i] = 0.;
    }
}


void 
solve_flow(double *restrict u, 
           double *restrict v, 
           double dx, 
           double dy,
           double *restrict p, 
           double *restrict rho,
           double *restrict nu,
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
        for ( i = 0; i < NX*NY; i++ ){
            un[i] = u[i];
            vn[i] = v[i];
        }

        solve_pressure_poisson(p, dx, dy, dt, u, v, rho);
        solve_stokes_momentum(u, v, un, vn, p, rho, nu, dt, dx, dy);

        apply_vel_boundary_conditions(u, v);
        printf("sc: %d\n", stepcount);
        // Check if in steady state
        udif = 0.;
        vdif = 0.;
        unt = 0.;
        vnt = 0.;

        #pragma omp parallel for
        for ( i = 0; i < NX*NY; i++ ){
            udif += fabs(fabs(u[i]) - fabs(un[i]));
            unt += fabs(un[i]);
            
            vdif += fabs(fabs(v[i]) - fabs(vn[i]));
            vnt += fabs(vn[i]);
        }
        if ( unt > 0. && vnt > 0. )
            diff = fabs( fabs(udif/unt) - fabs(vdif/vnt) ) / 2.;
        else
            diff = 1.0;

        stepcount++;
    } 
}
