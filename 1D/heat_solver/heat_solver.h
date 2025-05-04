/**
 * \defgroup heat_solver_1d 1D Heat Solver
 * \brief Functions related to 1D Heat Solver tests.
 *
 * \ingroup one_d
 */

/**
 * \ingroup heat_solver_1d
 * \file heat_solver.h
 * \brief Header for 1D heat equation solver.
 */

#ifndef HEAT_SOLVER_H
#define HEAT_SOLVER_H

#include <string>
#include "MDArray.h"

/**
 * @brief Compute Laplacian using central differences (Dirichlet BC).
 */
double laplacian(const MDArray<double>& u, size_t i, double dx);

/**
 * @brief Perform a single RK4 integration step using finite differences.
 */
void RK4_step(MDArray<double>& u, double dt, double dx, double alpha);

/**
 * @brief Save the solution array to a file for visualization.
 */
void saveSolution(const MDArray<double>& u, const std::string& filename, double dx, double time);

/**
 * @brief Perform a single RK4 step in Fourier space using the spectral method.
 */
void spectral_RK4_step(MDArray<double>& u, double dt, double alpha, double L);

/**
 * @brief Perform a single Backward Euler step in Fourier space using the spectral method.
 */
void spectral_BE_step(MDArray<double>& u, double dt, double alpha, double L);

#endif // HEAT_SOLVER_H
