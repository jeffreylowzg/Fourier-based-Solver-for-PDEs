/**
 * \defgroup convergence_1d 1D Convergence
 * \brief Functions related to 1D convergence tests.
 *
 * \ingroup one_d
 */

/**
 * \ingroup convergence_1d
 * \file heat_solver.h
 * \brief Header for 1D heat equation solver.
 */

/** heat_solver.h
 * 
 * Header file for 1D BE and RK4 code for testing convergence
 */

#ifndef HEAT_SOLVER_H
#define HEAT_SOLVER_H

#include <string>
#include "MDArray.h"

/**
 * @brief Saves the current solution array to a file for visualization.
 * 
 * @param u The solution array (MDArray<double>).
 * @param filename The name of the output file.
 * @param dx The spatial grid spacing.
 * @param time The current simulation time.
 * 
 * @details Writes the array values and their corresponding spatial positions to a file,
 * including the simulation time in the header.
 */
void saveSolution(const MDArray<double>& u, const std::string& filename, double dx, double time);

/**
 * @brief Advances the solution in time using the spectral method with RK4 integration.
 * 
 * @param u The solution array (MDArray<double>), updated in place.
 * @param dt The time step size.
 * @param alpha The thermal diffusivity coefficient.
 * @param L The length of the spatial domain.
 * @param step The number of time steps to perform.
 * 
 * @details Transforms the solution to Fourier space using FFT, applies RK4 integration 
 * to each Fourier mode, and transforms back to physical space.
 */
void spectral_RK4_steps(MDArray<double>& u, double dt, double alpha, double L, int step);

/**
 * @brief Advances the solution in time using the spectral method with Backward Euler integration.
 * 
 * @param u The solution array (MDArray<double>), updated in place.
 * @param dt The time step size.
 * @param alpha The thermal diffusivity coefficient.
 * @param L The length of the spatial domain.
 * @param step The number of time steps to perform.
 * 
 * @details Transforms the solution to Fourier space using FFT, applies Backward Euler 
 * integration to each Fourier mode, and transforms back to physical space.
 */
void spectral_BE_steps(MDArray<double>& u, double dt, double alpha, double L, int step);

#endif // HEAT_SOLVER_H
