/** heat_solver2d.h
 * 
 * header file 2D heat solving methods
 * includes fd RK4, spectral RK4, and spectral BE
 */

#ifndef HEAT_SOLVER2D_H
#define HEAT_SOLVER2D_H

#include <string>
#include "MDArray2D.h"

// Finite-difference methods (Dirichlet BC):
double laplacian2d(const MDArray2D<double>& u, size_t i, size_t j, double dx);
void RK4_step_2d(MDArray2D<double>& u, double dt, double dx, double alpha);
// RK4 with source:
void RK4_step_2d_source(MDArray2D<double>& u,
                        double dt, double dx,
                        double alpha,
                        double t, double L,
                        const std::string& source_type);

// Save solution for visualization:
void saveSolution2D(const MDArray2D<double>& u,
                    const std::string& filename,
                    double dx,
                    double time);

// Spectral methods (periodic BC over [0,L]x[0,L]):
//   pure diffusion:
void spectral_RK4_step_2d(MDArray2D<double>& u,
                          double dt, double alpha,
                          double L);
void spectral_BE_step_2d(MDArray2D<double>& u,
                         double dt, double alpha,
                         double L);
//   with source:
void spectral_RK4_step_2d_source(MDArray2D<double>& u,
                                 double dt, double alpha,
                                 double L, double t,
                                 const std::string& source_type);
void spectral_BE_step_2d_source(MDArray2D<double>& u,
                                double dt, double alpha,
                                double L, double t,
                                const std::string& source_type);

#endif // HEAT_SOLVER2D_H
