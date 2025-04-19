#include "heat_solver2d.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <fftw3.h>

//-------------------------------------------------------------------
// Spectral RK4 step using FFTW (Forward integration in Fourier space)
//-------------------------------------------------------------------
void spectral_RK4_step_2d(MDArray2D<double>& u, double dt, double alpha, double L) {
    int n = u.rows();
    int nc = n;  // n x n grid.
    // Allocate FFTW arrays.
    double* in = fftw_alloc_real(nc * nc);
    fftw_complex* out = fftw_alloc_complex(nc * (nc/2 + 1));
    
    // Copy u into in.
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            in[i*nc+j] = u(i, j);
        }
    }
    
    fftw_plan forward = fftw_plan_dft_r2c_2d(nc, nc, in, out, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);
    
    // Loop over Fourier modes.
    for (int i = 0; i < nc; i++) {
        int k_i = (i <= nc/2) ? i : i - nc;
        double kx = 2.0 * M_PI * k_i / L;
        for (int j = 0; j < nc/2 + 1; j++) {
            int k_j = j;
            double ky = 2.0 * M_PI * k_j / L;
            double ksq = kx * kx + ky * ky;
            double lambda = -alpha * ksq;
            
            int index = i * (nc/2 + 1) + j;
            double a = out[index][0];
            double b = out[index][1];
            
            // RK4 stages in Fourier space: dU/dt = lambda U.
            double k1_a = lambda * a;
            double k1_b = lambda * b;
            
            double a_k2 = a + 0.5 * dt * k1_a;
            double b_k2 = b + 0.5 * dt * k1_b;
            double k2_a = lambda * a_k2;
            double k2_b = lambda * b_k2;
            
            double a_k3 = a + 0.5 * dt * k2_a;
            double b_k3 = b + 0.5 * dt * k2_b;
            double k3_a = lambda * a_k3;
            double k3_b = lambda * b_k3;
            
            double a_k4 = a + dt * k3_a;
            double b_k4 = b + dt * k3_b;
            double k4_a = lambda * a_k4;
            double k4_b = lambda * b_k4;
            
            out[index][0] = a + dt / 6.0 * (k1_a + 2*k2_a + 2*k3_a + k4_a);
            out[index][1] = b + dt / 6.0 * (k1_b + 2*k2_b + 2*k3_b + k4_b);
        }
    }
    
    fftw_plan backward = fftw_plan_dft_c2r_2d(nc, nc, out, in, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);
    
    // Normalize and copy back into u.
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            u(i, j) = in[i*nc+j] / (nc * nc);
        }
    }
    fftw_free(in);
    fftw_free(out);
}

//-------------------------------------------------------------------
// Spectral Backward Euler step using FFTW (Backward integration in Fourier space)
//-------------------------------------------------------------------
void spectral_BE_step_2d(MDArray2D<double>& u, double dt, double alpha, double L) {
    int n = u.rows();
    int nc = n;
    double* in = fftw_alloc_real(nc * nc);
    fftw_complex* out = fftw_alloc_complex(nc * (nc/2 + 1));
    
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            in[i*nc+j] = u(i, j);
        }
    }
    
    fftw_plan forward = fftw_plan_dft_r2c_2d(nc, nc, in, out, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);
    
    // Update each Fourier mode with Backward Euler:
    // (1 - dt * lambda) * U^(n+1) = U^(n)  with lambda = -alpha*k^2, so 1 - dt*lambda = 1 + dt*alpha*k^2.
    for (int i = 0; i < nc; i++) {
        int k_i = (i <= nc/2) ? i : i - nc;
        double kx = 2.0 * M_PI * k_i / L;
        for (int j = 0; j < nc/2 + 1; j++) {
            int k_j = j;
            double ky = 2.0 * M_PI * k_j / L;
            double ksq = kx*kx + ky*ky;
            double factor = 1.0 / (1.0 + dt * alpha * ksq);
            int index = i * (nc/2 + 1) + j;
            out[index][0] *= factor;
            out[index][1] *= factor;
        }
    }
    
    fftw_plan backward = fftw_plan_dft_c2r_2d(nc, nc, out, in, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);
    
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            u(i, j) = in[i*nc+j] / (nc * nc);
        }
    }
    fftw_free(in);
    fftw_free(out);
}
