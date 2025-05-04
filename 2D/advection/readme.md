# 2D Advection (Convection) Solver
This is a 2D convection solver with full and pseudo spectral methods. This code can be run with backward Euler or RK4 and three different initial conditions.

Additionally, you can compare the full and pseudo RK4 results.

## Usage

You can run it with
```
  ./convection2d.exe <method> <initial_condition>
```
Where:

 `<method>` is the method to use. The options are:
 - spectral_rk4 // pseudo
 - spectral_be  // pseudo
 - full_spectral_rk4 // full
 - compare // runs both pseudo and full and saves results from both

`<initial_condition>` is the initial condition type:
 - sine
 - gaussian
 - disk
Both arguments are required

Example: `./convection2d.exe compare disk`
