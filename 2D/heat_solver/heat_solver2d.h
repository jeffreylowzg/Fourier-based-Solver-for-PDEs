/**
 * \defgroup heat_solver_2d 2D heat solver
 * \brief Functions related to heat solver
 *
 * \ingroup two_d
 */

/**
 * \ingroup heat_solver_2d
 * \file heat_solver2d.h
 * \brief Header for 2D heat solver. Includes fd RK4, spectral RK4, and spectral BE
 */

#ifndef HEAT_SOLVER2D_H
#define HEAT_SOLVER2D_H

#include <string>
#include "MDArray2D.h"

/// @brief Finite-difference methods (Dirichlet BC)
/// @param u input
/// @param i row
/// @param j col
/// @param dx grid spacing
/// @return the 2D laplacian
double laplacian2d(const MDArray2D<double>& u, size_t i, size_t j, double dx);

/// @brief rk4 step
/// @param u input
/// @param dt timestep
/// @param dx grid spacing
/// @param alpha diffusion coefficient
void RK4_step_2d(MDArray2D<double>& u, double dt, double dx, double alpha);

/// @brief RK4 with source
/// @param u input
/// @param dt timestep
/// @param dx grid spacing
/// @param alpha diffusion coefficient
/// @param t current time
/// @param L domain size
/// @param source_type the type of source
void RK4_step_2d_source(MDArray2D<double>& u,
                        double dt, double dx,
                        double alpha,
                        double t, double L,
                        const std::string& source_type);

/// @brief Save solution for visualization
/// @param u input
/// @param filename file to write to 
/// @param dx grid spacing
/// @param time time of step
void saveSolution2D(const MDArray2D<double>& u,
                    const std::string& filename,
                    double dx,
                    double time);

// Spectral methods (periodic BC over [0,L]x[0,L]):

/// @brief rk4 step
/// @param u input
/// @param dt timestep
/// @param alpha diffusion coefficient
/// @param L domain size
void spectral_RK4_step_2d(MDArray2D<double>& u,
                          double dt, double alpha,
                          double L);

/// @brief be step
/// @param u input
/// @param dt timestep
/// @param alpha diffusion coefficient
/// @param L domain size
void spectral_BE_step_2d(MDArray2D<double>& u,
                         double dt, double alpha,
                         double L);

/// @brief rk4 step with source
/// @param u input
/// @param dt timestep
/// @param alpha diffusion coefficient
/// @param L domain size
/// @param t time
/// @param source_type type of source
void spectral_RK4_step_2d_source(MDArray2D<double>& u,
                                 double dt, double alpha,
                                 double L, double t,
                                 const std::string& source_type);

/// @brief be step with source
/// @param u input
/// @param dt timestep
/// @param alpha diffusion coefficient
/// @param L domain size
/// @param t time
/// @param source_type type of source
void spectral_BE_step_2d_source(MDArray2D<double>& u,
                                double dt, double alpha,
                                double L, double t,
                                const std::string& source_type);

#endif // HEAT_SOLVER2D_H
