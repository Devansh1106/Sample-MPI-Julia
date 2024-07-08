push!(LOAD_PATH, "linhyp")
using Grid
using HypSys1D
using EqLinAdv

grid_size = 10              # number of cells
xmin, xmax = 0.0, 1.0       # domain
nvar = 3                    # number of variables
fprime = [1.0 2.0 3.0;
          0.0 -2.0 3.0;
          0.0 0.0 3.0]

final_time = 1.0
numflux = "rusanov"
boundary_condition = "Periodic"
initial_value(x) = [sin(2.0*pi*x),0.5*sin(2.0*pi*x),1.5*sin(2.0*pi*x)]
cfl = 0.9

equation = LinAdv(fprime)
problem = Problem((xmin, xmax), nvar, initial_value, boundary_condition, final_time)
param = create_parameters(cfl, grid_size)
scheme = Scheme(equation, numflux)
solve(equation, problem, scheme, param)
