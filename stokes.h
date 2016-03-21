#ifndef _stokes_h
#define _stokes_h

#include "environment.h"

void 
solve_pressure_poisson(double *p, 
                       double dx, 
                       double dy, 
                       double dt, 
                       double *u, 
                       double *v,
                       double *rho);

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
                      double dy);

void 
apply_vel_boundary_conditions(double *u, 
                              double *v);

void 
solve_flow(double *u, 
           double *v, 
           double dx, 
           double dy,
           double *p, 
           double *rho,
           double *nu,
           double dt);

#endif
