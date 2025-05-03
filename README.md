# 18847A Final Project: Spectral Solver for PDEs

Authors:
 - Jeffrey Low
 - Jacob Rohozen

This directory contains code for solving heat equations with pseudo-spectral and full-spectral methods in 1D and 2D. In 2D, there is code to step forward systems with advection and with source terms. We use Backward Euler and RK4 to step the system forward. Additionally, our presentation slides are included in the directory.

File structure:

1D  
 - convergence  
 - heat_solver  

2D
 - advection  
 - heat_solver  
 - reaction  

presentation_slides
 - .pdf
 - .pptx

## Building and Running

In the top-level directory, edit the Make.include settings to find your FFTW.  
Navigate to the solver you would like to try and run `make`. This will build the executable. `make run` will run the executable with the default settings.
If you would like to try different settings (integrators, initial conditions, sources) run the executable and select from the settings it prints out.

Each file has more information available at the top.

## Doxygen
Install doxygen

```sudo apt install doxygen```

In the root of this git run the following to generate the doxygen documentation

```doxygen Doxyfile```

Open the following file in html to view the generated documentation

```Fourier-based-Solver-for-PDEs/html/index.html```
