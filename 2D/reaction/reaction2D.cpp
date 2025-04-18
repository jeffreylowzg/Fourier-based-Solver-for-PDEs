#include "reaction2D.h"
#include <fftw3.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/stat.h>

//-------------------------------------------------------------------
// Helper: Apply Dirichlet boundary conditions (u = 0 on the boundary).
//-------------------------------------------------------------------
static void applyBoundaryConditions(MDArray2D<double>& u) {
    int n = static_cast<int>(u.rows());
    for (int i = 0; i < n; i++) {
        u(i, 0) = 0.0;
        u(i, n - 1) = 0.0;
    }
    for (int j = 0; j < n; j++) {
        u(0, j) = 0.0;
        u(n - 1, j) = 0.0;
    }
}

//-------------------------------------------------------------------
// Diffusion update using a spectral RK4 method (for the diffusion term only).
// The linear operator is L(k) = -D * k^2, with k^2 = kx^2 + ky^2.
//-------------------------------------------------------------------
static void diffusion_RK4_step(MDArray2D<double>& u, double dt, double D, double L) {
    int n = static_cast<int>(u.rows());
    int nc = n;  // Assume square grid.
    double* in = fftw_alloc_real(nc * nc);
    fftw_complex* out = fftw_alloc_complex(nc * (nc/2 + 1));

    // Copy u into the FFTW input array.
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            in[i * nc + j] = u(i, j);
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
            int index = i * (nc/2 + 1) + j;
            double ky = 2.0 * M_PI * j / L;
            double ksq = kx * kx + ky * ky;
            if (ksq == 0.0) // Skip the zero mode.
                continue;
            double a = out[index][0];
            double b = out[index][1];
            double lambda = -D * ksq;
            
            // RK4 stages (update is applied to both real and imaginary parts).
            double k1a = lambda * a;
            double k1b = lambda * b;
            
            double a2 = a + 0.5 * dt * k1a;
            double b2 = b + 0.5 * dt * k1b;
            double k2a = lambda * a2;
            double k2b = lambda * b2;
            
            double a3 = a + 0.5 * dt * k2a;
            double b3 = b + 0.5 * dt * k2b;
            double k3a = lambda * a3;
            double k3b = lambda * b3;
            
            double a4 = a + dt * k3a;
            double b4 = b + dt * k3b;
            double k4a = lambda * a4;
            double k4b = lambda * b4;
            
            out[index][0] = a + dt / 6.0 * (k1a + 2*k2a + 2*k3a + k4a);
            out[index][1] = b + dt / 6.0 * (k1b + 2*k2b + 2*k3b + k4b);
        }
    }
    
    fftw_plan backward = fftw_plan_dft_c2r_2d(nc, nc, out, in, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);
    
    double invN2 = 1.0 / (nc * nc);
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            u(i, j) = in[i * nc + j] * invN2;
        }
    }
    
    fftw_free(in);
    fftw_free(out);
    applyBoundaryConditions(u);
}

//-------------------------------------------------------------------
// Diffusion update using a spectral backward Euler (BE) method.
// In Fourier space the update is:
//    u_hat^(n+1) = u_hat^(n) / (1 + dt * D * k^2).
//-------------------------------------------------------------------
static void diffusion_BE_step(MDArray2D<double>& u, double dt, double D, double L) {
    int n = static_cast<int>(u.rows());
    int nc = n;
    double* in = fftw_alloc_real(nc * nc);
    fftw_complex* out = fftw_alloc_complex(nc * (nc/2 + 1));
    
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            in[i * nc + j] = u(i, j);
        }
    }
    
    fftw_plan forward = fftw_plan_dft_r2c_2d(nc, nc, in, out, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);
    
    for (int i = 0; i < nc; i++) {
        int k_i = (i <= nc/2) ? i : i - nc;
        double kx = 2.0 * M_PI * k_i / L;
        for (int j = 0; j < nc/2 + 1; j++) {
            int index = i * (nc/2 + 1) + j;
            double ky = 2.0 * M_PI * j / L;
            double ksq = kx * kx + ky * ky;
            if (ksq == 0.0)
                continue;
            double a = out[index][0];
            double b = out[index][1];
            double denom = 1.0 + dt * D * ksq;
            out[index][0] = a / denom;
            out[index][1] = b / denom;
        }
    }
    
    fftw_plan backward = fftw_plan_dft_c2r_2d(nc, nc, out, in, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);
    
    double invN2 = 1.0 / (nc * nc);
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            u(i, j) = in[i * nc + j] * invN2;
        }
    }
    
    fftw_free(in);
    fftw_free(out);
    applyBoundaryConditions(u);
}

//-------------------------------------------------------------------
// Reaction update using a pointwise RK4 method for f(u) = r u (1-u).
//-------------------------------------------------------------------
static void reaction_RK4_update(MDArray2D<double>& u, double dt, double r) {
    int n = static_cast<int>(u.rows());
    int m = static_cast<int>(u.cols());
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double u0 = u(i, j);
            double k1 = r * u0 * (1 - u0);
            double k2 = r * (u0 + 0.5 * dt * k1) * (1 - (u0 + 0.5 * dt * k1));
            double k3 = r * (u0 + 0.5 * dt * k2) * (1 - (u0 + 0.5 * dt * k2));
            double k4 = r * (u0 + dt * k3) * (1 - (u0 + dt * k3));
            u(i, j) = u0 + dt / 6.0 * (k1 + 2*k2 + 2*k3 + k4);
        }
    }
}

//-------------------------------------------------------------------
// Reaction update using a simple explicit Euler method.
//-------------------------------------------------------------------
static void reaction_Euler_update(MDArray2D<double>& u, double dt, double r) {
    int n = static_cast<int>(u.rows());
    int m = static_cast<int>(u.cols());
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double u0 = u(i, j);
            u(i, j) = u0 + dt * r * u0 * (1 - u0);
        }
    }
}

//-------------------------------------------------------------------
// Reaction–Diffusion step using Strang splitting (RK4 version).
// Performs: half diffusion (RK4) -> full reaction (RK4) -> half diffusion (RK4).
//-------------------------------------------------------------------
void reactiondiffusion_RK4_step_2d(MDArray2D<double>& u,
                                   double dt, double D, double r,
                                   double L) {
    diffusion_RK4_step(u, dt / 2.0, D, L);
    reaction_RK4_update(u, dt, r);
    diffusion_RK4_step(u, dt / 2.0, D, L);
}

//-------------------------------------------------------------------
// Reaction–Diffusion step using Strang splitting (BE version for diffusion).
// Performs: half diffusion (BE) -> full reaction (Euler) -> half diffusion (BE).
//-------------------------------------------------------------------
void reactiondiffusion_BE_step_2d(MDArray2D<double>& u,
                                  double dt, double D, double r,
                                  double L) {
    diffusion_BE_step(u, dt / 2.0, D, L);
    reaction_Euler_update(u, dt, r);
    diffusion_BE_step(u, dt / 2.0, D, L);
}

//-------------------------------------------------------------------
// Save the 2D solution to a file.
// The file format includes time, grid coordinates, and u values.
//-------------------------------------------------------------------
void saveSolution2D(const MDArray2D<double>& u,
                    const std::string& filename, double dx, double time) {
    std::ofstream file("data/" + filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open file data/" << filename << " for writing." << std::endl;
        return;
    }
    file << "# Time: " << time << "\n";
    size_t rows = u.rows();
    size_t cols = u.cols();
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double x = i * dx;
            double y = j * dx;
            file << x << " " << y << " " << u(i, j) << "\n";
        }
        file << "\n";
    }
    file.close();
}
