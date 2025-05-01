#ifndef HEAT_SOLVER_H
#define HEAT_SOLVER_H

#include <string>
#include "MDArray.h"

//-------------------------------------------------------------------
// Function declarations for helper routines (FD and spectral)
//-------------------------------------------------------------------
double laplacian(const MDArray<double>& u, size_t i, double dx);
void RK4_step(MDArray<double>& u, double dt, double dx, double alpha);
void saveSolution(const MDArray<double>& u, const std::string& filename, double dx, double time);

// Spectral methods (assumes periodic BC)
// Forward spectral solver using RK4 integration in Fourier space
void spectral_RK4_step(MDArray<double>& u, double dt, double alpha, double L);

// Backward spectral solver using Backward Euler integration in Fourier space
void spectral_BE_step(MDArray<double>& u, double dt, double alpha, double L);

#endif // HEAT_SOLVER_H
