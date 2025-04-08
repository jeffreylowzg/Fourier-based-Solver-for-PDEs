#ifndef SPECTRAL_ADVECTION_DIFFUSION_H
#define SPECTRAL_ADVECTION_DIFFUSION_H

#include <string>
#include "MDArray2D.h"

// Spectral RK4 solver for the 2D convection–diffusion (now "convection")
// equation: u_t + v·grad(u) = D ∇²u, with imposed Dirichlet boundaries.
void spectral_RK4_step_2d_convection(MDArray2D<double>& u,
                                     double dt, double D,
                                     double vx, double vy,
                                     double L);

// Spectral Backward Euler solver for the 2D convection–diffusion equation,
// with imposed Dirichlet boundaries.
void spectral_BE_step_2d_convection(MDArray2D<double>& u,
                                    double dt, double D,
                                    double vx, double vy,
                                    double L);

void saveSolution2D(const MDArray2D<double>& u, const std::string& filename, double dx, double time);


#endif // SPECTRAL_ADVECTION_DIFFUSION_H
