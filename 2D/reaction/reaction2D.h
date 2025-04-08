#ifndef REACTION_DIFFUSION_H
#define REACTION_DIFFUSION_H

#include <string>
#include "MDArray2D.h"

//-------------------------------------------------------------------
// Reaction–Diffusion solvers using Strang splitting.
// Two methods are provided:
//  1. reactiondiffusion_RK4_step_2d: Uses a spectral RK4 solver for diffusion
//     and an RK4 update for the nonlinear reaction term.
//  2. reactiondiffusion_BE_step_2d: Uses a spectral backward Euler solver for diffusion
//     and an explicit Euler update for the reaction term.
//
// The PDE solved is:
//    u_t = D ∇² u + r u (1 - u)
// with Dirichlet (zero) boundary conditions.
//-------------------------------------------------------------------
void reactiondiffusion_RK4_step_2d(MDArray2D<double>& u,
                                   double dt, double D, double r,
                                   double L);

void reactiondiffusion_BE_step_2d(MDArray2D<double>& u,
                                  double dt, double D, double r,
                                  double L);

// Save the 2D solution to a file (for visualization).
void saveSolution2D(const MDArray2D<double>& u,
                    const std::string& filename, double dx, double time);

#endif // REACTION_DIFFUSION_H
