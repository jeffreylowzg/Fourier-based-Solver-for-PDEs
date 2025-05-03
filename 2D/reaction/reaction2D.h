/**
 * \defgroup reaction_2d 2D reaction solver
 * \brief Functions related to 2D reaction
 *
 * \ingroup two_d
 */

/**
 * \ingroup reaction_2d
 * \file reaction2D.h
 * \brief Header for 2D reaction solver.
 */

#ifndef REACTION_DIFFUSION_H
#define REACTION_DIFFUSION_H

#include <string>
#include "MDArray2D.h"

//-------------------------------------------------------------------
// Reaction-Diffusion solvers using Strang splitting.
// Two methods are provided:
//  1. reactiondiffusion_RK4_step_2d: spectral RK4 for diffusion,
//     RK4 for reaction (no source).
//  2. reactiondiffusion_BE_step_2d: spectral backward Euler for diffusion,
//     explicit Euler for reaction (no source).
//
// The PDE solved is:
//    u_t = D grad^2 u + r u (1 - u)
// with homogeneous Dirichlet boundary conditions.
//-------------------------------------------------------------------

/// @brief Strang split: RK4 diffusion + RK4 reaction (no source)
/// @param u input
/// @param dt timestep
/// @param D diffusion coefficient
/// @param r reaction coefficient
/// @param L domain size
void reactiondiffusion_RK4_step_2d(MDArray2D<double>& u,
                                   double dt,
                                   double D,
                                   double r,
                                   double L);

/// @brief Strang split: BE diffusion + Euler reaction (no source)
/// @param u input
/// @param dt timestep
/// @param D diffusion coefficient
/// @param r reaction coefficient
/// @param L domain size
void reactiondiffusion_BE_step_2d(MDArray2D<double>& u,
                                  double dt,
                                  double D,
                                  double r,
                                  double L);

//-------------------------------------------------------------------
// Strang-splitting solvers with selectable source term S(x,y,t).
// Pass the current time and a source-type string.
//-------------------------------------------------------------------

/// @brief Strang split + source: RK4 diffusion + RK4 reaction+source
/// @param u input
/// @param dt timestep
/// @param D diffusion coefficient
/// @param r reaction coefficient
/// @param L domain size
/// @param current_time the current simulation time
/// @param source_type type of source
void reactiondiffusion_RK4_step_2d_source(MDArray2D<double>& u,
                                          double dt,
                                          double D,
                                          double r,
                                          double L,
                                          double current_time,
                                          const std::string& source_type);

/// @brief Strang split + source: BE diffusion + RK4 reaction+source
/// @param u input
/// @param dt timestep
/// @param D diffusion coefficient
/// @param r reaction coefficient
/// @param L domain size
/// @param current_time the current simulation time
/// @param source_type type of source
void reactiondiffusion_BE_step_2d_source(MDArray2D<double>& u,
                                         double dt,
                                         double D,
                                         double r,
                                         double L,
                                         double current_time,
                                         const std::string& source_type);

/// @brief RK4 reaction diffusion with model and source
/// @param u reactant 1
/// @param v reactant 2
/// @param dt timestep
/// @param D diffusion coefficient
/// @param r reaction coefficient
/// @param L domain size
/// @param t time
/// @param source_type type of source
/// @param model_type type of model (logistic or grayscott)
void reactiondiffusion_RK4_step_2d_combined(MDArray2D<double>& u,
                                            MDArray2D<double>& v,
                                            double dt, double D, double r,
                                            double L, double t,
                                            const std::string& source_type,
                                            const std::string& model_type);

/// @brief backward Euler reaction diffusion with model and source
/// @param u reactant 1
/// @param v reactant 2
/// @param dt timestep
/// @param D diffusion coefficient
/// @param r reaction coefficient
/// @param L domain size
/// @param t time
/// @param source_type type of source
/// @param model_type type of model (logistic or grayscott)
void reactiondiffusion_BE_step_2d_combined(MDArray2D<double>& u,
                                           MDArray2D<double>& v,
                                           double dt, double D, double r,
                                           double L, double t,
                                           const std::string& source_type,
                                           const std::string& model_type);


/// @brief Save the 2D solution to a file (for visualization).
/// @param u heat grid
/// @param filename output file
/// @param dx grid spacing
/// @param time time in simulation
void saveSolution2D(const MDArray2D<double>& u,
                    const std::string& filename,
                    double dx,
                    double time);

#endif // REACTION_DIFFUSION_H
