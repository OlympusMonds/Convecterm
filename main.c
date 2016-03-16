#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Domain size
#define XMIN 0.
#define XMAX 5.
#define YMIN 0.
#define YMAX 2.5

// Resolution
#define NX 101
#define NY 51

#define GRAVITY 0.0003

// Characters for visualisation
#define COLOR_ESC    '\033'
#define COLOR_RESET  "\033[0m"
#define BOLD         "\033[1m"
#define WIPE_SCREEN  "\033[2J"


void visualise(double (*arr)[NX], double min, double max, int clear){
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
           if ( arr[j][i] < r1 )
               field[count+3] = '4';  // BLUE
           else if (arr[j][i] < r2 && arr[j][i] >= r1 )
               field[count+3] = '6';  // CYAN
           else if (arr[j][i] < r3 && arr[j][i] >= r2 )
               field[count+3] = '2';  // GREEN
           else if (arr[j][i] < r4 && arr[j][i] >= r3 )
               field[count+3] = '7';  // WHITE
           else if (arr[j][i] < r5 && arr[j][i] >= r4 )
               field[count+3] = '3';  // YELLOW
           else if (arr[j][i] < r6 && arr[j][i] >= r5 )
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


void solve_pressure_poisson(double (*p)[NX], double dx, double dy, 
                            double dt, double (*u)[NX], double (*v)[NX],
                            double (*rho)[NX]){
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


void solve_stokes_momentum(double (*u)[NX], double (*v)[NX],
                           double (*un)[NX], double (*vn)[NX],
                           double (*p)[NX], double (*rho)[NX],
                           double (*nu)[NX],
                           double dt, double dx, double dy){
    int i, j;

    // Pre-compute
    double dtodx2 = dt / (dx * dx);
    double dtody2 = dt / (dy * dy);
   
    for ( j = 1; j < NY-1; j++ ){
        for ( i = 1; i < NX-1; i++){
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


void apply_thermal_boundary_conditions(double (*t)[NX]){
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


void apply_vel_boundary_conditions(double (*u)[NX], double (*v)[NX]){
    int i, j;
    
    for ( j = 0; j < NY; j++ ){
        // Vel Left wall, freeslip
        u[j][0] = 0.;
        v[j][0] = v[j][1];

        // Right wall, freeslip
        u[j][NX-1] = 0.;
        v[j][NX-1] = v[j][NX-2];
    }

    for ( i = 6; i < NX-6; i++ ){
        // Bottom wall, freeslip
        u[0][i] = u[1][i];
        v[0][i] = 0.;

        // Top wall, freeslip
        u[NY-1][i] = u[NY-2][i];
        v[NY-1][i] = 0.;
    }
}


void solve_advection_diffusion(double (*t)[NX], double (*u)[NX], double (*v)[NX],
                               double dx, double dy,
                               double (*rho)[NX], double dt,
                               double cp, double (*k)[NX],
                               double H){
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


void solve_flow(double (*u)[NX], double (*v)[NX], 
                double dx, double dy,
                double (*p)[NX], double (*rho)[NX],
                double (*nu)[NX],
                double dt){

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
            diff = ( fabs(udif/unt) - fabs(vdif/vnt) ) / 2.;
        else
            diff = 1.0;

        stepcount++;
    } 
}


void update_nu(double (*nu)[NX], double (*t)[NX]){
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


void update_rho(double (*rho)[NX], double (*t)[NX]){
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


void update_k(double (*k)[NX], double (*t)[NX]){
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


void calc_dt(double (*u)[NX], double (*v)[NX], double (*rho)[NX],
             double* dx, double* dy, 
             double* dt, double* k, 
             double* cp){
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
           uv = fabs(u[j][i]);
           vv = fabs(v[j][i]);
           max_velx = uv > max_velx ? uv : max_velx;
           max_vely = vv > max_vely ? vv : max_vely;
           min_rho = rho[j][i] < min_rho ? rho[j][i] : min_rho;
       }
    }
    max_vel = max_velx > max_vely ? max_velx : max_vely;

    min_seperation = *dx < *dy ? *dx : *dy;
    if ( max_vel != 0 )
        vel_dt = 1e-3 * min_seperation / (max_vel*2);

    nrg_dt =  (min_seperation * min_seperation) / (*k / min_rho * *cp);

    *dt = vel_dt < nrg_dt ? vel_dt : nrg_dt;
}


int main () {

    double current_time;

    double dx = (XMAX - XMIN) / (NX - 1);
    double dy = (YMAX - YMIN) / (NY - 1);

    //double x[NX];
    //double y[NY];
    double u[NY][NX];    // vel in x
    double v[NY][NX];    // vel in y
    double p[NY][NX];    // pressure
    double rho[NY][NX];  // density
    double nu[NY][NX];   // viscosity
    double t[NY][NX];    // temperature
    double k[NY][NX];    // conductivity

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
          u[j][i] = 0.; 
          v[j][i] = 0.; 
          p[j][i] = ref_rho * dy * GRAVITY * (NY - j);  // Approximate lithostatic pressure
          rho[j][i] = ref_rho; 
          nu[j][i] = 1.;
          // Make the temp field unstable 
          if ( j < 4 && i < (int)(NX/2) )
              t[j][i] = 1000.;
          else if ( j > NY-5 && i > (int)(NX/2) )
              t[j][i] = 0.;
          else
              t[j][i] = 500.;
          k[j][i] = 100.;
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
