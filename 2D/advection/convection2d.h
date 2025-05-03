/**
 * \defgroup convection_2d 2D convection solver
 * \brief Functions related to 2d convection
 *
 * \ingroup two_d
 */

/**
 * \ingroup convection_2d
 * \file convection2d.h
 * \brief Header for 2D convection solver.
 */

#ifndef SPECTRAL_ADVECTION_DIFFUSION_H
#define SPECTRAL_ADVECTION_DIFFUSION_H

#include <string>
#include "MDArray2D.h"

/// @brief Spectral RK4 step for the 2D convection-diffusion equation
/// @param u input
/// @param dt timestep
/// @param D diffusion constant
/// @param vx velocity field x
/// @param vy velocity field y
/// @param L domain length
void spectral_RK4_step_2d_convection(MDArray2D<double>& u,
                                     double dt, double D,
                                     double vx, double vy,
                                     double L);

/// @brief Spectral Backward Euler step for the 2D convection-diffusion equation
/// @param u input
/// @param dt timestep
/// @param D diffusion constant
/// @param vx velocity field x
/// @param vy velocity field y
/// @param L domain length
void spectral_BE_step_2d_convection(MDArray2D<double>& u,
                                    double dt, double D,
                                    double vx, double vy,
                                    double L);

/// @brief Save the 2D solution to a file
/// @param u input
/// @param filename file to save to
/// @param dx grid step size
/// @param time current time
void full_spectral_RK4_step_2d_convection(MDArray2D<double>& u,
                                          double dt, double D,
                                          double vx, double vy,
                                          double L);

/// @brief Save the 2D solution to a file
/// @param u input
/// @param filename file to save to
/// @param dx grid step size
/// @param time current time
void saveSolution2D(const MDArray2D<double>& u, const std::string& filename, double dx, double time);


#endif // SPECTRAL_ADVECTION_DIFFUSION_H
