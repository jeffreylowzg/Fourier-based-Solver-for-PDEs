/** heat_solver.h
 * 
 * header file for 1D BE and RK4 code for testing convergence
 */

#ifndef HEAT_SOLVER_H
#define HEAT_SOLVER_H

#include <string>
#include "MDArray.h"

//-------------------------------------------------------------------
// Function declarations for helper routines (FD and spectral)
//-------------------------------------------------------------------
void saveSolution(const MDArray<double>& u, const std::string& filename, double dx, double time);

// Spectral methods (assumes periodic BC)
// Forward spectral solver using RK4 integration in Fourier space
void spectral_RK4_steps(MDArray<double>& u, double dt, double alpha, double L, int step);

// Backward spectral solver using Backward Euler integration in Fourier space
void spectral_BE_steps(MDArray<double>& u, double dt, double alpha, double L, int step);

#endif // HEAT_SOLVER_H
