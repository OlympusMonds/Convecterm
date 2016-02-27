#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define XMIN 0.
#define XMAX 2.
#define YMIN 0.
#define YMAX 2.

#define NX 51
#define NY 51

#define NIT 100
#define NT 100000
#define MAX_TIME 2000.
#define GRAVITY 0.0001

#define COLOR_RESET  "\033[0m"
#define BOLD         "\033[1m"
#define BLACK_TEXT   "\033[30;1m"
#define RED_TEXT     "\033[31;1m"
#define GREEN_TEXT   "\033[32;1m"
#define YELLOW_TEXT  "\033[33;1m"
#define BLUE_TEXT    "\033[34;1m"
#define MAGENTA_TEXT "\033[35;1m"
#define CYAN_TEXT    "\033[36;1m"
#define WHITE_TEXT   "\033[37;1m"
#define WIPE_SCREEN  "\033[2J"

void visualise(double arr[][NX], double min, double max){
    int i,j;

    double low = ((max - min) * 0.1429) + min;
    double mlow = ((max - min) * 0.28571) + min;
    double mlow2 = ((max - min) * 0.42857) + min;
    double mid = ((max - min) * 0.57128) + min;
    double mhigh = ((max - min) * 0.7142857) + min;
    double high = ((max - min) * 0.85714) + min;
    
    printf(WIPE_SCREEN);
    for ( j = NY-1; j >= 0; j-- ){
       for ( i = 0; i < NX; i++ ){
           if ( arr[j][i] < low )
               printf(BLUE_TEXT);
           else if (arr[j][i] < mlow && arr[j][i] >= low )
               printf(CYAN_TEXT);
           else if (arr[j][i] < mlow2 && arr[j][i] >= mlow )
               printf(GREEN_TEXT);
           else if (arr[j][i] < mid && arr[j][i] >= mlow2 )
               printf(WHITE_TEXT);
           else if (arr[j][i] < mhigh && arr[j][i] >= mid )
               printf(YELLOW_TEXT);
           else if (arr[j][i] < high && arr[j][i] >= mhigh )
               printf(RED_TEXT);
           else
               printf(MAGENTA_TEXT);
           printf("\u2588\u2588");
       }
       printf("\n");
    }
    printf(COLOR_RESET);
}

void solve_pressure_poisson(double p[][NX], double* dx, double* dy, 
                            double *dt, double u[][NX], double v[][NX],
                            double rho[][NX]){
    double b[NY][NX];
    double pn[NY][NX];

    int i, j, n;
    double dx2 = *dx * *dx;
    double dy2 = *dy * *dy;

    for ( j = 1; j < NY-1; j++ ){
        for ( i = 1; i < NX-1; i++){
            b[j][i] = (( rho[j][i] * dx2 * dy2 ) / ( 2 * (dx2 + dy2))) * 
                      ( 
                          1 / *dt * 
                          (
                           (( u[j][i+1] - u[j][i-1] ) / ( 2 * *dx )) +
                           (( v[j+1][i] - v[j-1][i] ) / ( 2 * *dy ))
                          ) -
                          pow( (( u[j][i+1] - u[j][i-1] ) / ( 2* *dx)), 2) -
                          2 *
                          ((( u[j+1][i] - u[j-1][i] ) / ( 2 * *dy )) *
                           (( v[j][i+1] - v[j][i-1] ) / ( 2 * *dx ))) -
                          pow( (( v[j+1][i] - v[j-1][i] ) / ( 2 * *dy )), 2)
                      );
        }
    } 

    for ( n = 0; n < NIT; n++ ){
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
    }



}


void solve_stokes_momentum(double u[][NX], double v[][NX],
                           double un[][NX], double vn[][NX],
                           double p[][NX], double rho[][NX],
                           double nu[][NX],
                           double *dt, double *dx, double *dy){
    int i, j;

    for ( j = 1; j < NY-1; j++ ){
        for ( i = 1; i < NX-1; i++){
            u[j][i] = un[j][i] - ( *dt / (rho[j][i] * 2. * *dx) ) * (p[j][i+1] - p[j][i-1]) +
                      nu[j][i] * (
                                  ((*dt/ (*dx * *dx)) * (un[j][i+1] - 2*un[j][i] + un[j][i-1])) +
                                  ((*dt/ (*dy * *dy)) * (un[j+1][i] - 2*un[j][i] + un[j-1][i]))
                                 );
            
            v[j][i] = (vn[j][i] - ( *dt / (rho[j][i] * 2. * *dy) ) * (p[j+1][i] - p[j-1][i]) +
                      nu[j][i] * (
                                  ((*dt/ (*dx * *dx)) * (vn[j][i+1] - 2*vn[j][i] + vn[j][i-1])) +
                                  ((*dt/ (*dy * *dy)) * (vn[j+1][i] - 2*vn[j][i] + vn[j-1][i]))
                                 )) - (GRAVITY * rho[j][i]);
        }
    }
}


void apply_thermal_boundary_conditions(double t[][NX]){
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


void apply_vel_boundary_conditions(double u[][NX], double v[][NX]){
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
        // Bottom wall, imposed shear
        //u[0][i] = -1e-1;
        u[0][i] = u[1][i];
        v[0][i] = 0.;

        // Top wall, imposed shear
        //u[NY-1][i] = 1e-1;
        u[NY-1][i] = u[NY-2][i];
        v[NY-1][i] = 0.;
    }
}


void solve_advection_diffusion(double t[][NX], double u[][NX], double v[][NX],
                               double* dx, double* dy,
                               double rho[][NX], double* dt,
                               double* cp, double* k,
                               double* H){
    int i,j; 
    double tn[NY][NX];
    double kx;
    double ky;

    for ( j = 0; j < NY; j++ ){
       for ( i = 0; i < NX; i++ ){
           tn[j][i] = t[j][i];
       }
    }
    
    for ( j = 1; j < NY-1; j++ ){
        for ( i = 1; i < NX-1; i++){
           kx = *k * (tn[j][i+1] - 2.*tn[j][i] + tn[j][i-1]) / (*dx * *dx);
           ky = *k * (tn[j+1][i] - 2.*tn[j][i] + tn[j-1][i]) / (*dy * *dy);

           t[j][i] = tn[j][i] + *dt * ((*H + kx + ky)/(rho[j][i] * *cp) \
                     - (u[j][i] * ( (tn[j][i+1] - tn[j][i-1]) / (2 * *dx) )) \
                     - (v[j][i] * ( (tn[j+1][i] - tn[j-1][i]) / (2 * *dy) )) );
       }
    }

    apply_thermal_boundary_conditions(t);
}


void solve_flow(double u[][NX], double v[][NX], 
                double* dx, double* dy,
                double p[][NX], double rho[][NX],
                double nu[][NX],
                double* dt){

    double un[NY][NX];
    double vn[NY][NX];
    double diff = 1000.;
    double udif, vdif, unt, vnt;
    int stepcount = 0;
    int i, j;

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
               vdif += fabs(fabs(v[j][i]) - fabs(vn[j][i]));
               unt += fabs(un[j][i]);
               vnt += fabs(vn[j][i]);
           }
        }
        if ( unt != 0. && vnt != 0. )
            diff = ( fabs(udif/unt) - fabs(vdif/vnt) ) / 2.;
        else
            diff = 1.0;

        stepcount++;
        //printf("res: %e, sc: %d\n",diff, stepcount);
    } 
}


void update_nu(double nu[][NX], double t[][NX]){
    int i, j;

    double ref_nu = 1.;
    double ref_temp = 500.;
    double theta = 0.1;

    for ( j = 0; j < NY; j++ ){
       for ( i = 0; i < NX; i++ ){
           nu[j][i] = ref_nu * exp(-theta * (t[j][i] - ref_temp)/ref_temp);
       }
    }
}


void update_rho(double rho[][NX], double t[][NX]){
    int i, j;

    double ref_rho = 100.;
    double ref_temp = 500.;
    double thermal_expansivity = 0.001;

    for ( j = 0; j < NY; j++ ){
       for ( i = 0; i < NX; i++ ){
           rho[j][i] = ref_rho * (1 - (thermal_expansivity * (t[j][i] - ref_temp)));
       }
    }
}


void calc_dt(double u[][NX], double v[][NX], double rho[][NX],
             double* dx, double* dy, 
             double* dt, double* k, 
             double* cp){
    int i,j;
    double vel_dt = 1e6;;
    double nrg_dt;
    double min_seperation;
    double max_vel;
    double min_rho;
    double uv, vv;

    max_vel = 0.;
    min_rho = 1e6;

    for ( j = 0; j < NY; j++ ){
       for ( i = 0; i < NX; i++ ){
           uv = fabs(u[j][i]);
           vv = fabs(v[j][i]);
           max_vel = uv > max_vel ? uv : max_vel;
           max_vel = vv > max_vel ? vv : max_vel;  // TODO: split this into 2 vars, for potential compiler opt
           min_rho = rho[j][i] < min_rho ? rho[j][i] : min_rho;
       }
    }

    min_seperation = *dx < *dy ? *dx : *dy;
    if ( max_vel != 0 )
        vel_dt = 1e-3 * min_seperation / (max_vel*2);

    nrg_dt =  (min_seperation * min_seperation) / (*k / min_rho * *cp);

    //printf("maxvel: %f, min_rho: %f, min_sep: %f, vel_dt: %e, nrg_dt: %e", max_vel, min_rho, min_seperation, vel_dt, nrg_dt);

    //*dt = vel_dt < nrg_dt ? vel_dt : nrg_dt;
    *dt = nrg_dt;
    //printf(", dt = %e\n", *dt);
}


int main () {

    double dx = (XMAX - XMIN) / (NX - 1);
    double dy = (YMAX - YMIN) / (NY - 1);
    

    double x[NX];
    double y[NY];
    double u[NY][NX];
    double v[NY][NX];
    double p[NY][NX];
    double b[NY][NX];
    double rho[NY][NX];
    double nu[NY][NX];
    double t[NY][NX];

    double cp = 100.;
    double k = 100.;
    double H = 0.;
    double dt = 5e-5;

    int i, j;  // Init some counters

    for ( i = 0; i < NX; i++ ) {
        x[i] = XMIN + i * dx;
    }
    for ( i = 0; i < NY; i++ ) {
        y[i] = YMIN + i * dy;
    }

    // Initial conditions
    for ( j = 0; j < NY; j++ ){
       for ( i = 0; i < NX; i++ ){
          u[j][i] = 0.; 
          v[j][i] = 0.; 
          p[j][i] = 0.; 
          b[j][i] = 0.;
          rho[j][i] = 100.; 
          nu[j][i] = 1.; 
          if ( j < 4 && i < (int)(NX/2) ) {
              t[j][i] = 1000.;
              v[j][i] = 0.1; 
          }
          else if ( j > NY-4 && i > (int)(NX/2) ){
              t[j][i] = 0.;
              v[j][i] = -0.1; 
          }
          else
              t[j][i] = 500.;
       }
    }

    // Init the boundaries too
    apply_vel_boundary_conditions(u, v);
    apply_thermal_boundary_conditions(t);



    int timestep = 0;
    double current_time = 0;
    
    // Main loop
    //while ( timestep < NT && current_time < MAX_TIME ) {
    while ( 1 ) {
        if ( timestep % 1000 == 0 ) {
            visualise(t, 0., 1000.);
            //visualise(u, -1e-2, 1e-2);
            printf("Timestep: %04d, Current time: %.4f, dt: %e\n", timestep, current_time, dt);
            printf("t: %3.3f, %3.3f, nu: %3.3f, %3.3f, rho: %3.3f, %3.3f\n", t[1][20],
                                                                             t[40][20],
                                                                             nu[1][20],
                                                                             nu[40][20],
                                                                             rho[1][20],
                                                                             rho[40][20]);
            printf("\n");
        }

        solve_advection_diffusion(t, u, v, &dx, &dy, rho, &dt, &cp, &k, &H);
        update_rho(rho, t);
        update_nu(nu, t);
        solve_flow(u, v, &dx, &dy, p, rho, nu, &dt);

        current_time += dt;
        timestep++;
        //calc_dt(u, v, rho, &dx, &dy, &dt, &k, &cp);
    }























    return 0;
}
