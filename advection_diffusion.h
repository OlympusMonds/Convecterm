#ifndef _advection_diffusion_h
#define _advection_diffusion_h

#include "environment.h"

void 
apply_thermal_boundary_conditions(double (*t)[NX]);

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
                          double H);

void 
update_nu(double (*nu)[NX], 
          double (*t)[NX]);


void 
update_rho(double (*rho)[NX], 
           double (*t)[NX]);


void 
update_k(double (*k)[NX], 
         double (*t)[NX]);


#endif
