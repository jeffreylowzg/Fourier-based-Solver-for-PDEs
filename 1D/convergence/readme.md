# 1D Convergence

This file tests convergence by stepping the simulation forward in time and comparing the result to the analytical solution.

The initial condition is a sine wave with period 1 on the domain \[0, 1\].  
The analytical solution is a decaying exponential.

The convergence test runs multiple times; you can vary the grid spacing and/or the timestep.

## Usage

```
./convergence_test <test> [options]
```

where:

 `<test>` is either `be` or `rk4`,  

`[options]` are `vary_n` and/or `vary_dt`.

- `vary_n` → halves the grid spacing after each test (finer grid)  
- `vary_dt` → halves the timestep after each test (smaller timestep)
