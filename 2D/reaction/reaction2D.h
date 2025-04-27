#ifndef REACTION_DIFFUSION_H
#define REACTION_DIFFUSION_H

#include <string>
#include "MDArray2D.h"

//-------------------------------------------------------------------
// Reaction–Diffusion solvers using Strang splitting.
// Two methods are provided:
//  1. reactiondiffusion_RK4_step_2d: spectral RK4 for diffusion,
//     RK4 for reaction (no source).
//  2. reactiondiffusion_BE_step_2d: spectral backward Euler for diffusion,
//     explicit Euler for reaction (no source).
//
// The PDE solved is:
//    u_t = D ∇²u + r u (1 - u)
// with homogeneous Dirichlet boundary conditions.
//-------------------------------------------------------------------
void reactiondiffusion_RK4_step_2d(MDArray2D<double>& u,
                                   double dt,
                                   double D,
                                   double r,
                                   double L);

void reactiondiffusion_BE_step_2d(MDArray2D<double>& u,
                                  double dt,
                                  double D,
                                  double r,
                                  double L);

//-------------------------------------------------------------------
// Strang‐splitting solvers with selectable source term S(x,y,t).
// Pass the current time and a source‐type string.
//-------------------------------------------------------------------
void reactiondiffusion_RK4_step_2d_source(MDArray2D<double>& u,
                                          double dt,
                                          double D,
                                          double r,
                                          double L,
                                          double current_time,
                                          const std::string& source_type);

void reactiondiffusion_BE_step_2d_source(MDArray2D<double>& u,
                                         double dt,
                                         double D,
                                         double r,
                                         double L,
                                         double current_time,
                                         const std::string& source_type);

void reactiondiffusion_RK4_step_2d_combined(MDArray2D<double>& u,
                                            MDArray2D<double>& v,
                                            double dt, double D, double r,
                                            double L, double t,
                                            const std::string& source_type,
                                            const std::string& model_type);

void reactiondiffusion_BE_step_2d_combined(MDArray2D<double>& u,
                                           MDArray2D<double>& v,
                                           double dt, double D, double r,
                                           double L, double t,
                                           const std::string& source_type,
                                           const std::string& model_type);


//-------------------------------------------------------------------
// Save the 2D solution to a file (for visualization).
//-------------------------------------------------------------------
void saveSolution2D(const MDArray2D<double>& u,
                    const std::string& filename,
                    double dx,
                    double time);

#endif // REACTION_DIFFUSION_H
