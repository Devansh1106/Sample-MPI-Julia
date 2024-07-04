# Solver for Riemann problem of Linear Advection
# User chooses Ul, Ur
custom_path = "/home/devansh/Sample-MPI-Julia/linhyp"
push!(LOAD_PATH,custom_path)
using Grid
using FV
using EqLinAdv
using LinearAlgebra
using Plots
# plotly(size = (750, 565)) # use plotly for interactivity

grid_size = 50             # number of cells
xmin, xmax = 0.0, 1.0      # domain   
nvar = 3                   # number of variables
fprime(U,x,eq)  = [1.0 2.0 3.0;
                   0.0 -2.0 3.0;
                   0.0 0.0 3.0] # The matrix given by F'(U)
                # The last argument is absolutely dummy to solve
                # Euler with the same code. Should we make `eq` an
                # optional argument?
final_time = 1.0

numflux   = "rusanov"

# initial condition
# Ul, Ur       = [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]
# initial_value(x) = (x <= 0) ? Ul : Ur
# boundary_value(x) = 0.0 # Dummy
boundary_condition = "Periodic"
initial_value(x) = [sin(2.0*pi*x), 0.5*sin(2.0*pi*x), 1.5*sin(2.0*pi*x)]

# save_time_interval = final_time
# skip_plotting = false
cfl = 0.0 # Currently not used
Ccfl = 0.9

# -----------------------------------
# plotters = get_plot_funcs(skip_plotting)
equation = get_equation(fprime)
problem = Problem((xmin,xmax), nvar, initial_value, boundary_value,
                  boundary_condition, final_time)
param = Parameters(grid_size, cfl, Ccfl, save_time_interval)
scheme = Scheme(equation, numflux)
@time plt_data = solve(equation, problem, scheme, param, plotters)
if plt_data !== nothing
   p, anim = plt_data
   savefig(p, "final_soln.png")
   gif(anim, "soln.gif", fps = 5) # would have been better in the solve function
                        # here because of VS Code
   plot(p, legend=true) # final solution
end

# TODO - compare with characteristic pictures in Ch 3