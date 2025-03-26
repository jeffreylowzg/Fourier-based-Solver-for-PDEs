// heat_solver.cpp
#include "heat_solver.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <sys/stat.h>
#include <fftw3.h>

//-------------------------------------------------------------------
// Compute Laplacian using central differences (Dirichlet BC)
//-------------------------------------------------------------------
double laplacian(const MDArray<double>& u, size_t i, double dx) {
    if(i == 0 || i == u.size()-1)
        return 0.0;
    return (u[i+1] - 2*u[i] + u[i-1]) / (dx * dx);
}

//-------------------------------------------------------------------
// RK4 (Explicit) integration step (finite-difference)
//-------------------------------------------------------------------
void RK4_step(MDArray<double>& u, double dt, double dx, double alpha) {
    size_t n = u.size();
    MDArray<double> k1(n), k2(n), k3(n), k4(n), temp(n);
    
    for(size_t i = 0; i < n; i++) 
        k1[i] = alpha * laplacian(u, i, dx);

    for(size_t i = 0; i < n; i++) 
        temp[i] = u[i] + 0.5 * dt * k1[i];
    for(size_t i = 0; i < n; i++) 
        k2[i] = alpha * laplacian(temp, i, dx);

    for(size_t i = 0; i < n; i++) 
        temp[i] = u[i] + 0.5 * dt * k2[i];
    for(size_t i = 0; i < n; i++) 
        k3[i] = alpha * laplacian(temp, i, dx);

    for(size_t i = 0; i < n; i++) 
        temp[i] = u[i] + dt * k3[i];
    for(size_t i = 0; i < n; i++) 
        k4[i] = alpha * laplacian(temp, i, dx);

    for(size_t i = 0; i < n; i++)
        u[i] += dt / 6.0 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}

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
// Spectral step using FFTW (for periodic BC)
//-------------------------------------------------------------------
void spectral_step(MDArray<double>& u, double dt, double alpha, double L) {
    int n = u.size();
    // Allocate memory for FFTW arrays
    double* in = fftw_alloc_real(n);
    fftw_complex* out = fftw_alloc_complex(n/2 + 1);

    // Copy the current solution into the input array
    for (int i = 0; i < n; i++) {
        in[i] = u[i];
    }

    // Create and execute forward FFT (real-to-complex)
    fftw_plan forward = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);

    // Update Fourier coefficients exactly for the heat equation:
    // The wave number for mode k is: k_val = 2*pi*k/L
    // and û(k, t+dt) = û(k, t) * exp(-alpha * k_val^2 * dt)
    for (int k = 0; k < n/2 + 1; k++) {
        double k_val = (2.0 * M_PI * k) / L;
        double factor = exp(-alpha * k_val * k_val * dt);
        out[k][0] *= factor; // real part
        out[k][1] *= factor; // imaginary part
    }

    // Create and execute inverse FFT (complex-to-real)
    fftw_plan backward = fftw_plan_dft_c2r_1d(n, out, in, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);

    // Normalize the result (FFTW does not scale the inverse transform)
    for (int i = 0; i < n; i++) {
        u[i] = in[i] / n;
    }

    // Free FFTW allocated memory
    fftw_free(in);
    fftw_free(out);
}
