#ifndef BENCH
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "stokes.h"
#include "advection_diffusion.h"


void 
visualise(double *arr, 
          double min, 
          double max, 
          int clear)
{
    /* Visualise a 2D field, coloured with 7 colours between min 
     * and max. If clear is 0, then clear the screen before visualising,
     */

    int i,j;
    unsigned int count;

    // 7 colours in the terminal, split into bins
    double r1 = ((max - min) * 0.1429) + min;
    double r2 = ((max - min) * 0.28571) + min;
    double r3 = ((max - min) * 0.42857) + min;
    double r4 = ((max - min) * 0.57128) + min;
    double r5 = ((max - min) * 0.71429) + min;
    double r6 = ((max - min) * 0.85714) + min;

    // Store characters in field before printfing.
    char field[(NY * (9*NX + 1)) + 1];

    if ( clear != 0 )
       printf(WIPE_SCREEN);

    printf(BOLD);

    count = 0;
    for ( j = NY-1; j >= 0; j-- ) {  // Origin is bottom-left
       for ( i = 0; i < NX; i++ ) {
           field[count] = COLOR_ESC;
           field[count+1] = '[';
           field[count+2] = '3';
           if ( arr[NX*j + i] < r1 )
               field[count+3] = '4';  // BLUE
           else if (arr[NX*j + i] < r2 && arr[NX*j + i] >= r1 )
               field[count+3] = '6';  // CYAN
           else if (arr[NX*j + i] < r3 && arr[NX*j + i] >= r2 )
               field[count+3] = '2';  // GREEN
           else if (arr[NX*j + i] < r4 && arr[NX*j + i] >= r3 )
               field[count+3] = '7';  // WHITE
           else if (arr[NX*j + i] < r5 && arr[NX*j + i] >= r4 )
               field[count+3] = '3';  // YELLOW
           else if (arr[NX*j + i] < r6 && arr[NX*j + i] >= r5 )
               field[count+3] = '1';  // RED
           else
               field[count+3] = '5';  // MAGENTA
           field[count+4] = ';';
           field[count+5] = '1';
           field[count+6] = 'm';
           field[count+7] = '#';
           field[count+8] = '#'; 
           count += 9;
       }
       field[count] = '\n';
       count += 1;
    }
    field[count] = '\0';
    printf("%s", field);
    printf(COLOR_RESET);
}


void 
calc_dt(double *u, 
        double *v, 
        double *rho,
        double* dx, 
        double* dy, 
        double* dt, 
        double* k, 
        double* cp)
{
    int i,j;
    double vel_dt = 1e6;
    double nrg_dt;
    double min_seperation;
    double max_vel, max_velx, max_vely;
    double min_rho;
    double uv, vv;

    max_velx = 0.;
    max_vely = 0.;
    min_rho = 1e6;

    for ( j = 0; j < NY; j++ ){
       for ( i = 0; i < NX; i++ ){
           uv = fabs(u[NX*j + i]);
           vv = fabs(v[NX*j + i]);
           max_velx = uv > max_velx ? uv : max_velx;
           max_vely = vv > max_vely ? vv : max_vely;
           min_rho = rho[NX*j + i] < min_rho ? rho[NX*j + i] : min_rho;
       }
    }
    max_vel = max_velx > max_vely ? max_velx : max_vely;

    min_seperation = *dx < *dy ? *dx : *dy;
    if ( max_vel != 0 )
        vel_dt = 1e-3 * min_seperation / (max_vel*2);

    nrg_dt =  (min_seperation * min_seperation) / (*k / min_rho * *cp);

    *dt = vel_dt < nrg_dt ? vel_dt : nrg_dt;
}


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
    

    // Don't actually use these guys
    //for ( i = 0; i < NX; i++ ) {
    //    x[i] = XMIN + i * dx;
    //}
    //for ( i = 0; i < NY; i++ ) {
    //    y[i] = YMIN + i * dy;
    //}

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
    while ( 1 ) {
        if ( timestep % 1000 == 0 ) {
            visualise(t, 0., 1000., 1);
            printf("Timestep: %06d, Current time: %.4f, dt: %e\n", timestep, current_time, dt);
            printf("\n");
        }

        // Solve thermal stuff first, so the flow equations have some meat to start with
        solve_advection_diffusion(t, u, v, dx, dy, rho, dt, cp, k, H);
        update_rho(rho, t);
        update_nu(nu, t);
        update_k(k, t);
        solve_flow(u, v, dx, dy, p, rho, nu, dt);

        current_time += dt;
        timestep++;
        //calc_dt(u, v, rho, &dx, &dy, &dt, &k, &cp);  // TODO: model is very sensitive to ANY dt change. Fix this.
    }


    return 0;
}
#endif
