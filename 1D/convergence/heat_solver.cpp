/** heat_solver.cpp
 * 
 * This file contains the 1D BE and RK4 code for testing convergence
 */

#include "heat_solver.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <sys/stat.h>
#include <fftw3.h>

//-------------------------------------------------------------------
// Save the solution to a file for visualization at given timestep
//-------------------------------------------------------------------
void saveSolution(const MDArray<double>& u, const std::string& filename, double dx, double time) {
    std::ofstream file("data/" + filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open file data/" << filename << " for writing." << std::endl;
        return;
    }
    file << "# Time: " << time << std::endl;
    for(size_t i = 0; i < u.size(); i++){
        double x = i * dx;
        file << x << " " << u[i] << "\n";
    }
    file.close();
}

//-------------------------------------------------------------------
// Spectral steps using RK4 integration (forward solver)
//-------------------------------------------------------------------
void spectral_RK4_steps(MDArray<double>& u, double dt, double alpha, double L, int steps) {
    int n = u.size();
    // Allocate FFTW arrays.
    double* in = fftw_alloc_real(n);
    fftw_complex* out = fftw_alloc_complex(n/2 + 1);

    // Forward FFT: copy the solution into the FFTW input array.
    for (int i = 0; i < n; i++) {
        in[i] = u[i];
    }
    fftw_plan forward = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);

    // For each Fourier mode, apply RK4 to solve: dŪ/dt = -alpha * k^2 * Ū.
    for (int step = 0; step < steps; step++) {
      for (int k = 0; k < n/2 + 1; k++) {
          double k_val = (2.0 * M_PI * k) / L;
          double lambda = -alpha * k_val * k_val;

          // Represent the Fourier coefficient as U = a + i*b.
          double a = out[k][0];
          double b = out[k][1];

          // RK4 stages:
          // k1 = lambda * U
          double k1_a = lambda * a;
          double k1_b = lambda * b;

          // k2 = lambda * (U + 0.5*dt*k1)
          double a_k2 = a + 0.5 * dt * k1_a;
          double b_k2 = b + 0.5 * dt * k1_b;
          double k2_a = lambda * a_k2;
          double k2_b = lambda * b_k2;

          // k3 = lambda * (U + 0.5*dt*k2)
          double a_k3 = a + 0.5 * dt * k2_a;
          double b_k3 = b + 0.5 * dt * k2_b;
          double k3_a = lambda * a_k3;
          double k3_b = lambda * b_k3;

          // k4 = lambda * (U + dt*k3)
          double a_k4 = a + dt * k3_a;
          double b_k4 = b + dt * k3_b;
          double k4_a = lambda * a_k4;
          double k4_b = lambda * b_k4;

          // Update the Fourier coefficient.
          out[k][0] = a + dt/6.0 * (k1_a + 2*k2_a + 2*k3_a + k4_a);
          out[k][1] = b + dt/6.0 * (k1_b + 2*k2_b + 2*k3_b + k4_b);
      }
    }

    // Inverse FFT to return to physical space.
    fftw_plan backward = fftw_plan_dft_c2r_1d(n, out, in, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);

    // Normalize and update u.
    for (int i = 0; i < n; i++) {
        u[i] = in[i] / n;
    }

    fftw_free(in);
    fftw_free(out);
}

//-------------------------------------------------------------------
// Spectral step using Backward Euler integration (backward solver)
//-------------------------------------------------------------------
void spectral_BE_steps(MDArray<double>& u, double dt, double alpha, double L, int steps) {
    int n = u.size();
    double* in = fftw_alloc_real(n);
    fftw_complex* out = fftw_alloc_complex(n/2 + 1);

    // Forward FFT.
    for (int i = 0; i < n; i++) {
        in[i] = u[i];
    }
    fftw_plan forward = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);

    for (int step = 0; step < steps; step++) {
      // Update each Fourier mode using backward Euler:
      // (1 + dt*alpha*k^2) * Ū^(n+1) = Ū^(n)  =>  Ū^(n+1) = Ū^(n) / (1 + dt*alpha*k^2)
      for (int k = 0; k < n/2 + 1; k++) {
          double k_val = (2.0 * M_PI * k) / L;
          double factor = 1.0 / (1.0 + alpha * k_val * k_val * dt);
          out[k][0] *= factor;
          out[k][1] *= factor;
      }
    }

    // Inverse FFT.
    fftw_plan backward = fftw_plan_dft_c2r_1d(n, out, in, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);

    // Normalize and update u.
    for (int i = 0; i < n; i++) {
        u[i] = in[i] / n;
    }

    fftw_free(in);
    fftw_free(out);
}
