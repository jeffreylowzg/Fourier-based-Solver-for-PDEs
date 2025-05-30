/** reaction2D.cpp
 * 
 * This file contains the implementation for 2D reaction using pseudo-spectral
 * methods.
 * 
 * Several model types, initial conditions, and source terms are supported:
 * 
 * model types:
 *  - logistic
 *  - grayscott
 * 
 * initial conditions:
 *  - gaussian
 *  - gaussian16
 *  - zero
 * 
 * source terms:
 *  - none
 *  - gaussian_pulse
 *  - standing_wave
 *  - traveling_gaussian
 *  - pulsed_ring
 * 
 * The solver works by using Strang splitting:
 *  1. diffusion step with dt/2
 *  2. reaction step with dt
 *  3. diffusion step with dt/2
 */

#include "reaction2D.h"
#include <algorithm>      // for std::max, std::min
#include <fftw3.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <string>

/// @brief Source term definitions (standing_wave clamped >= 0)
/// @param x x position
/// @param y y position
/// @param t time
/// @param L domain size
/// @param source_type type of source
/// @return source value
double source_term(double x, double y, double t,
                   double L, const std::string& source_type) {
    double cx = 0.5 * L;
    double cy = 0.5 * L;

    if (source_type == "gaussian_pulse") {
        double A     = 1.0;
        double sigma = 0.1 * L;
        double tau   = 1.0;
        double dx2   = (x - cx)*(x - cx);
        double dy2   = (y - cy)*(y - cy);
        double spatial  = std::exp(-(dx2 + dy2)/(2.0*sigma*sigma));
        double temporal = 1.0 - std::exp(-t/tau);
        return A * spatial * temporal;
    }
    else if (source_type == "standing_wave") {
        double A      = 1.0;
        double period = 5.0;
        double omega  = 2.0 * M_PI / period;
        int mx = 4, my = 4;
        double s = A
                 * std::sin(mx*M_PI*x/L)
                 * std::sin(my*M_PI*y/L)
                 * std::cos(omega*t);
        // clamp negative lobes to zero
        return std::max(0.0, s);
    }
    else if (source_type == "traveling_gaussian") {
        double A     = 1.0;
        double v     = 0.2 * L;
        double sigma = 0.05 * L;
        double x0    = 0.1 * L;
        double xc    = x0 + v * t;
        double dx2   = (x - xc)*(x - xc);
        double dy2   = (y - cy)*(y - cy);
        return A * std::exp(-(dx2 + dy2)/(2.0*sigma*sigma));
    }
    else if (source_type == "pulsed_ring") {
        double r0      = std::sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy));
        double r1      = 0.2 * L;
        double r2      = 0.3 * L;
        double t0      = 1.0;
        double dt_p    = 0.5;
        if (r0 >= r1 && r0 <= r2 && t >= t0 && t <= t0 + dt_p)
            return 1.0;
        else
            return 0.0;
    }
    // none or unrecognized type
    return 0.0;
}

/// @brief dirichlet boundary conditions (u = 0 on boundary).
/// @param u input
static void applyBoundaryConditions(MDArray2D<double>& u) {
    int n = static_cast<int>(u.rows());
    for (int i = 0; i < n; i++) {
        u(i, 0)     = 0.0;
        u(i, n - 1) = 0.0;
        u(0, i)     = 0.0;
        u(n - 1, i) = 0.0;
    }
}

/// @brief Diffusion update using spectral RK4 (diffusion only).
/// @param u input
/// @param dt timestep
/// @param D diffusion coefficient
/// @param L domain size
static void diffusion_RK4_step(MDArray2D<double>& u,
                               double dt, double D, double L) {
    int nc = static_cast<int>(u.rows());
    double* in  = fftw_alloc_real(nc * nc);
    fftw_complex* out = fftw_alloc_complex(nc * (nc/2 + 1));

    for (int i = 0; i < nc; i++)
      for (int j = 0; j < nc; j++)
        in[i*nc + j] = u(i,j);

    auto plan_f = fftw_plan_dft_r2c_2d(nc, nc, in, out, FFTW_ESTIMATE);
    fftw_execute(plan_f);
    fftw_destroy_plan(plan_f);

    // RK4 in Fourier space
    for (int i = 0; i < nc; i++) {
        int k_i = (i <= nc/2) ? i : i - nc;
        double kx = 2.0 * M_PI * k_i / L;
        for (int j = 0; j < nc/2 + 1; j++) {
            int idx = i*(nc/2+1)+j;
            double ky  = 2.0 * M_PI * j / L;
            double ksq = kx*kx + ky*ky;
            if (ksq == 0.0) continue;
            double lambda = -D * ksq;
            double a = out[idx][0], b = out[idx][1];
            double k1a = lambda*a, k1b = lambda*b;
            double a2 = a + 0.5*dt*k1a, b2 = b + 0.5*dt*k1b;
            double k2a = lambda*a2, k2b = lambda*b2;
            double a3 = a + 0.5*dt*k2a, b3 = b + 0.5*dt*k2b;
            double k3a = lambda*a3, k3b = lambda*b3;
            double a4 = a + dt*k3a,     b4 = b + dt*k3b;
            double k4a = lambda*a4, k4b = lambda*b4;
            out[idx][0] = a + dt/6.0*(k1a + 2*k2a + 2*k3a + k4a);
            out[idx][1] = b + dt/6.0*(k1b + 2*k2b + 2*k3b + k4b);
        }
    }

    auto plan_b = fftw_plan_dft_c2r_2d(nc, nc, out, in, FFTW_ESTIMATE);
    fftw_execute(plan_b);
    fftw_destroy_plan(plan_b);

    double inv = 1.0/(nc*nc);
    for (int i = 0; i < nc; i++)
      for (int j = 0; j < nc; j++)
        u(i,j) = in[i*nc + j] * inv;

    fftw_free(in);
    fftw_free(out);
    applyBoundaryConditions(u);
}

/// @brief Diffusion update using spectral backward Euler (diffusion only).
/// @param u input
/// @param dt timestep
/// @param D diffusion coefficient
/// @param L domain size
static void diffusion_BE_step(MDArray2D<double>& u,
                              double dt, double D, double L) {
    int nc = static_cast<int>(u.rows());
    double* in  = fftw_alloc_real(nc * nc);
    fftw_complex* out = fftw_alloc_complex(nc * (nc/2 + 1));

    for (int i = 0; i < nc; i++)
      for (int j = 0; j < nc; j++)
        in[i*nc + j] = u(i,j);

    auto plan_f = fftw_plan_dft_r2c_2d(nc, nc, in, out, FFTW_ESTIMATE);
    fftw_execute(plan_f);
    fftw_destroy_plan(plan_f);

    for (int i = 0; i < nc; i++) {
        int k_i = (i <= nc/2) ? i : i - nc;
        double kx = 2.0 * M_PI * k_i / L;
        for (int j = 0; j < nc/2 + 1; j++) {
            int idx = i*(nc/2+1)+j;
            double ky = 2.0 * M_PI * j / L;
            double ksq = kx*kx + ky*ky;
            if (ksq == 0.0) continue;
            double factor = 1.0/(1.0 + dt*D*ksq);
            out[idx][0] *= factor;
            out[idx][1] *= factor;
        }
    }

    auto plan_b = fftw_plan_dft_c2r_2d(nc, nc, out, in, FFTW_ESTIMATE);
    fftw_execute(plan_b);
    fftw_destroy_plan(plan_b);

    double inv = 1.0/(nc*nc);
    for (int i = 0; i < nc; i++)
      for (int j = 0; j < nc; j++)
        u(i,j) = in[i*nc + j] * inv;

    fftw_free(in);
    fftw_free(out);
    applyBoundaryConditions(u);
}

/// @brief Reaction update using RK4 (no source)
/// @param u input
/// @param dt timestep
/// @param r reaction coefficient
static void reaction_RK4_update(MDArray2D<double>& u,
                                double dt, double r) {
    int n = static_cast<int>(u.rows());
    int m = static_cast<int>(u.cols());
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++) {
        double u0 = u(i,j);
        double k1 = r*u0*(1-u0);
        double k2 = r*(u0+0.5*dt*k1)*(1-(u0+0.5*dt*k1));
        double k3 = r*(u0+0.5*dt*k2)*(1-(u0+0.5*dt*k2));
        double k4 = r*(u0+    dt*k3)*(1-(u0+    dt*k3));
        u(i,j) = u0 + dt/6.0*(k1 + 2*k2 + 2*k3 + k4);
      }
}

/// @brief Reaction+source update using RK4, then clamp u in range [0,1]
/// @param u input
/// @param dt timestep
/// @param r reaction coefficient
/// @param t time
/// @param L domain size
/// @param src type of source
static void reaction_RK4_update_with_source(MDArray2D<double>& u,
                                            double dt, double r,
                                            double t, double L,
                                            const std::string& src) {
    int n = static_cast<int>(u.rows());
    int m = static_cast<int>(u.cols());
    double dx = L / n;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++) {
        double x = i*dx, y = j*dx, u0 = u(i,j);
        double s0 = source_term(x,y,t,      L, src);
        double k1 = r*u0*(1-u0) + s0;
        double th = t + 0.5*dt;
        double u2 = u0 + 0.5*dt*k1;
        double s1 = source_term(x,y,th,     L, src);
        double k2 = r*u2*(1-u2) + s1;
        double u3 = u0 + 0.5*dt*k2;
        double k3 = r*u3*(1-u3) + s1;
        double uf = u0 + dt*k3;
        double s2 = source_term(x,y,t+dt,   L, src);
        double k4 = r*uf*(1-uf) + s2;
        u(i,j) = u0 + dt/6.0*(k1 + 2*k2 + 2*k3 + k4);
      }
    // clamp into [0,1]
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < m; ++j)
        u(i,j) = std::max(0.0, std::min(1.0, u(i,j)));
}

/// @brief Strang split: RK4 diffusion + RK4 reaction (no source)
/// @param u input
/// @param dt timestep
/// @param D diffusion coefficient
/// @param r reaction coefficient
/// @param L domain size
void reactiondiffusion_RK4_step_2d(MDArray2D<double>& u,
                                   double dt, double D, double r,
                                   double L) {
    diffusion_RK4_step(u, dt/2.0, D, L);
    reaction_RK4_update(u, dt, r);
    diffusion_RK4_step(u, dt/2.0, D, L);
}

/// @brief Strang split: BE diffusion + Euler reaction (no source)
/// @param u input
/// @param dt timestep
/// @param D diffusion coefficient
/// @param r reaction coefficient
/// @param L domain size
void reactiondiffusion_BE_step_2d(MDArray2D<double>& u,
                                  double dt, double D, double r,
                                  double L) {
    diffusion_BE_step(u, dt/2.0, D, L);
    // no-source Euler reaction
    {
        int n = static_cast<int>(u.rows());
        int m = static_cast<int>(u.cols());
        for (int i = 0; i < n; i++)
          for (int j = 0; j < m; j++)
            u(i,j) += dt * r * u(i,j) * (1 - u(i,j));
    }
    diffusion_BE_step(u, dt/2.0, D, L);
}

/// @brief Strang split + source: RK4 diffusion + RK4 reaction+source
/// @param u input
/// @param dt timestep
/// @param D diffusion coefficient
/// @param r reaction coefficient
/// @param L domain size
/// @param current_time the current simulation time
/// @param source_type type of source
void reactiondiffusion_RK4_step_2d_source(MDArray2D<double>& u,
                                          double dt, double D, double r,
                                          double L, double t,
                                          const std::string& src) {
    diffusion_RK4_step(u, dt/2.0, D, L);
    reaction_RK4_update_with_source(u, dt, r, t, L, src);
    diffusion_RK4_step(u, dt/2.0, D, L);
}

/// @brief Strang split + source: BE diffusion + RK4 reaction+source
/// @param u input
/// @param dt timestep
/// @param D diffusion coefficient
/// @param r reaction coefficient
/// @param L domain size
/// @param current_time the current simulation time
/// @param source_type type of source
void reactiondiffusion_BE_step_2d_source(MDArray2D<double>& u,
                                         double dt, double D, double r,
                                         double L, double t,
                                         const std::string& src) {
    diffusion_BE_step(u, dt/2.0, D, L);
    reaction_Euler_update_with_source(u, dt, r, t, L, src);
    diffusion_BE_step(u, dt/2.0, D, L);
}

/// @brief Save the 2D solution to a file (for visualization).
/// @param u heat grid
/// @param filename output file
/// @param dx grid spacing
/// @param time time in simulation
void saveSolution2D(const MDArray2D<double>& u,
                    const std::string& filename,
                    double dx, double time) {
    std::ofstream file("data/" + filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open file data/" << filename << " for writing.\n";
        return;
    }
    file << "# Time: " << time << "\n";
    size_t rows = u.rows(), cols = u.cols();
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double x = i * dx, y = j * dx;
            file << x << " " << y << " " << u(i,j) << "\n";
        }
        file << "\n";
    }
    file.close();
}

/// @brief backward Euler reaction diffusion with model and source
/// @param u reactant 1
/// @param v reactant 2
/// @param dt timestep
/// @param D diffusion coefficient
/// @param r reaction coefficient
/// @param L domain size
/// @param t time
/// @param source_type type of source
/// @param model_type type of model (logistic or grayscott)
void reactiondiffusion_BE_step_2d_combined(MDArray2D<double>& u,
                                           MDArray2D<double>& v,
                                           double dt, double D, double r,
                                           double L, double t,
                                           const std::string& source_type,
                                           const std::string& model_type) {
    if (model_type == "logistic") {
        diffusion_BE_step(u, dt/2.0, D, L);
        reaction_Euler_update_with_source(u, dt, r, t, L, source_type);
        diffusion_BE_step(u, dt/2.0, D, L);
    }
    else if (model_type == "grayscott") {
        MDArray2D<double> u0 = u;
        MDArray2D<double> v0 = v;

        diffusion_BE_step(u0, dt/2.0, D, L);
        diffusion_BE_step(v0, dt/2.0, D, L);

        const double F = 0.04;
        const double k = 0.06;

        for (int i = 0; i < (int)u.rows(); ++i)
            for (int j = 0; j < (int)u.cols(); ++j) {
                double uu = u0(i,j), vv = v0(i,j);
                double uvv = uu * vv * vv;
                u(i,j) = uu + dt * (-uvv + F * (1.0 - uu));
                v(i,j) = vv + dt * ( uvv - (F + k) * vv);
            }

        diffusion_BE_step(u, dt/2.0, D, L);
        diffusion_BE_step(v, dt/2.0, D, L);
    }
    else {
        std::cerr << "Unknown model_type in BE step: " << model_type << "\n";
    }
}

/// @brief RK4 reaction diffusion with model and source
/// @param u reactant 1
/// @param v reactant 2
/// @param dt timestep
/// @param D diffusion coefficient
/// @param r reaction coefficient
/// @param L domain size
/// @param t time
/// @param source_type type of source
/// @param model_type type of model (logistic or grayscott)
void reactiondiffusion_RK4_step_2d_combined(MDArray2D<double>& u,
                                            MDArray2D<double>& v,
                                            double dt, double D, double r,
                                            double L, double t,
                                            const std::string& source_type,
                                            const std::string& model_type) {
    if (model_type == "logistic") {
        diffusion_RK4_step(u, dt/2.0, D, L);
        if (source_type != "none")
            reaction_RK4_update_with_source(u, dt, r, t, L, source_type);
        else
            reaction_RK4_update(u, dt, r);
        diffusion_RK4_step(u, dt/2.0, D, L);
    }
    else if (model_type == "grayscott") {
        MDArray2D<double> u0 = u;
        MDArray2D<double> v0 = v;

        diffusion_RK4_step(u0, dt/2.0, D, L);
        diffusion_RK4_step(v0, dt/2.0, D, L);

        const double F = 0.04;
        const double k = 0.06;

        for (int i = 0; i < (int)u.rows(); ++i)
            for (int j = 0; j < (int)u.cols(); ++j) {
                double uu = u0(i,j), vv = v0(i,j);
                double uvv = uu * vv * vv;
                u(i,j) = uu + dt * (-uvv + F * (1.0 - uu));
                v(i,j) = vv + dt * ( uvv - (F + k) * vv);
            }

        diffusion_RK4_step(u, dt/2.0, D, L);
        diffusion_RK4_step(v, dt/2.0, D, L);
    }
    else {
        std::cerr << "Unknown model_type in RK4 step: " << model_type << "\n";
    }
}
