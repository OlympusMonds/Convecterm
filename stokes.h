#ifndef _stokes_h
#define _stokes_h

#include "environment.h"

void 
solve_pressure_poisson(double (*p)[NX], 
                       double dx, 
                       double dy, 
                       double dt, 
                       double (*u)[NX], 
                       double (*v)[NX],
                       double (*rho)[NX]);

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
                      double dy);

void 
apply_vel_boundary_conditions(double (*u)[NX], 
                              double (*v)[NX]);

void 
solve_flow(double (*u)[NX], 
           double (*v)[NX], 
           double dx, 
           double dy,
           double (*p)[NX], 
           double (*rho)[NX],
           double (*nu)[NX],
           double dt);

#endif
