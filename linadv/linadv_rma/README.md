## Usage
To run the code in `linadv1D_rma_act.jl` and in `linadv1d_rma_pass.jl`:
```bash
cd linadv/linadv_rma
mpiexecjl -np <processes> julia <file-name.jl> <N> <Tt> <cfl> <a>
./rm.sh
```
- `<N>`: Total number of grid points

- `<Tf>`: Final time

- `<cfl>`: Courant–Friedrichs–Lewy (cfl) condition

- `<a>`: advection velocity

## Visualization
To visualize the results, run the following command after generating the solution files and running the `rm.sh` script file:

```bash
julia plot_rma.jl
```

