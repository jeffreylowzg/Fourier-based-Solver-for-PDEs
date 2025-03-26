#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cassert>
#include <string>
#include <sstream>
#include <sys/stat.h> 

//-------------------------------------------------------------------
// MDarray class for 1D arrays
//-------------------------------------------------------------------
template<typename T>
class MDArray {
private:
    size_t size_;
    std::vector<T> data;
public:
    MDArray(size_t n) : size_(n), data(n, T()) {}
    T& operator[](size_t i) { 
        assert(i < size_);
        return data[i]; 
    }
    const T& operator[](size_t i) const { 
        assert(i < size_);
        return data[i]; 
    }
    size_t size() const { return size_; }
};

//-------------------------------------------------------------------
// Compute Laplacian using central differences (Dirichlet BC)
//-------------------------------------------------------------------
double laplacian(const MDArray<double>& u, size_t i, double dx) {
    if(i == 0 || i == u.size()-1)
        return 0.0;
    return (u[i+1] - 2*u[i] + u[i-1]) / (dx * dx);
}

//-------------------------------------------------------------------
// RK4 (Explicit) integration step
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
// Main: Solves heat equation, saves multiple solutions over time
//-------------------------------------------------------------------
int main() {
    const double L = 1.0;
    const size_t n = 101;
    const double dx = L / (n - 1);
    const double alpha = 0.01;
    const double dt = 0.0001;
    const size_t steps = 5000;
    const size_t save_interval = 100; // Save every 100 steps

    mkdir("data", 0777);

    MDArray<double> u(n);

    // New initial condition: Gaussian pulse centered at L/2
    const double x0 = L / 2.0;   // Center of Gaussian
    const double sigma = 0.05;   // Width of Gaussian

    for(size_t i = 0; i < n; i++) {
        double x = i * dx;
        u[i] = exp(-( (x - x0)*(x - x0) ) / (2 * sigma * sigma));
    }

    // Initial condition save
    saveSolution(u, "solution_0.dat", dx, 0.0);

    // Time integration loop
    for(size_t step = 1; step <= steps; step++) {
        RK4_step(u, dt, dx, alpha);

        // Save at intervals
        if(step % save_interval == 0 || step == steps) {
            double current_time = step * dt;
            std::ostringstream fname;
            fname << "solution_" << step << ".dat";
            saveSolution(u, fname.str(), dx, current_time);
            std::cout << "Saved " << fname.str() << " at time t=" << current_time << std::endl;
        }
    }

    std::cout << "Simulation complete." << std::endl;
    return 0;
}
