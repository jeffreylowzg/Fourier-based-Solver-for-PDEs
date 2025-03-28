#ifndef HEAT_SOLVER2D_H
#define HEAT_SOLVER2D_H

#include <string>
#include "MDArray2D.h"

// Finite-difference methods:
double laplacian2d(const MDArray2D<double>& u, size_t i, size_t j, double dx);
void RK4_step_2d(MDArray2D<double>& u, double dt, double dx, double alpha);
void saveSolution2D(const MDArray2D<double>& u, const std::string& filename, double dx, double time);

// Spectral methods (assumes periodic boundary conditions over a square domain of side L):
// RK4 integration in Fourier space using FFTW.
void spectral_RK4_step_2d(MDArray2D<double>& u, double dt, double alpha, double L);
// Backward Euler integration in Fourier space using FFTW.
void spectral_BE_step_2d(MDArray2D<double>& u, double dt, double alpha, double L);

#endif // HEAT_SOLVER2D_H
