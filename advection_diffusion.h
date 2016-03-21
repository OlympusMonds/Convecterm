#ifndef _advection_diffusion_h
#define _advection_diffusion_h

#include "environment.h"

void 
apply_thermal_boundary_conditions(double *t);

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
                          double H);

void 
update_nu(double *nu, 
          double *t);


void 
update_rho(double *rho, 
           double *t);


void 
update_k(double *k, 
         double *t);


#endif
