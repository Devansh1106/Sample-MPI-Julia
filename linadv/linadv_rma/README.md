## Usage
To run the code in `linadv1D_rma_act.jl` and in `linadv1d_rma_pass.jl`:
```
cd linadv/linadv_rma
mpiexecjl -np <processes> julia <file-name.jl> <N> <Tt> <cfl> <a>
```
- `<N>`: Total number of grid points

- `<Tf>`: Final time

- `<cfl>`: Courant–Friedrichs–Lewy (cfl) condition

- `<a>`: advection velocity

## Visualization
To visualize the results, run the following command after generating the solution files:

```julia
julia plot.jl
```

