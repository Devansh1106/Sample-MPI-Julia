# solving linear advection equation 2D:-
# u_t + a*u_x + b*u_y= 0 with periodic bc taking a = b = 1
# initial conditions: 
# u(xy) = sin(2πx) * sin(2πy)
# Computational domain [0,1] × [0,1]

using DelimitedFiles
using Plots
using TimerOutputs
const to = TimerOutput()

# To generate a initial solution through initial condition
function initial_u!(param, x, y, u)
    u = @. sin(2.0 * π * x) * sin(2.0 * π * y)
end

# Lex-Wendroff method 
function update_lw!(u, unew, sigma_x, sigma_y)
    unew[1,1] = u[1,1] - 0.5 * sigma_x * (u[2,1] - u[end-1,1])
                       - 0.5 * sigma_y * (u[1,2] - u[1, end-1])
                       + 0.5 * sigma_x^2 * (u[end-1,1] - 2.0 * u[1,1] + u[2,1])
                       + 0.25 * sigma_x * sigma_y * (u[2,2] - u[end-1,2] - u[2,end-1] + u[end-1,end-1])
                       + 0.5 * sigma_y^2 * (u[1,end-1] - 2.0 * u[1,1] * u[1,2])
    @views unew[2:end-1,2:end-1] = u[2:end-1,2:end-1] - 0.5 * sigma_x * (u[3:end,2:end-1] - u[1:end-2,2:end-1])
                                                      - 0.5 * sigma_y * (u[2:end-1,3:end] - u[2:end-1,1:end-2])
                                                      + 0.5 * sigma_x^2 * (u[1:end-2,2:end-1] - 2.0*u[2:end-1,2:end-1] + u[3:end,2:end-1])
                                                      + 0.25 * sigma_x * sigma_y * (u[3:end,3:end] - u[1:end-2,3:end] - u[3:end,1:end] + u[1:end-2,1:end-2])
                                                      + 0.5 * sigma_y^2 * (u[2:end-1,1:end-2] - 2.0*u[2:end-1,2:end-1] + u[2:end-1,3:end])
    unew[end,end] = u[end,end] - 0.5 * sigma_x * (u[2,end] - u[end-1,end])
                               - 0.5 * sigma_y * (u[end,2] - u[end,end-1])
                               + 0.5 * sigma_x^2 * (u[end-1,end] - 2.0 * u[end,end] +u[2,end])
                               + 0.25 * sigma_x * sigma_y * (u[2,2] - u[end-1,2] - u[2,end-1] + u[end-1,end-1])
                               + 0.5 * sigma_y^2 * (u[end,end-1] - 2.0 * u[end,end] + u[end,2])
end

# Exact solution calculation
function exact_solution!(param, x, y, exact_sol)
    _a = param.a
    _Tf = param.Tf
    exact_sol = @. sin(2.0 * π * (x-(_a * _Tf))) * sin(2.0 * π * (y-(_a * _Tf)))
end

# Error calculation
function error_cal(param, exact_sol, u)
    error = 0.0
    error = sum(abs.(exact_sol - u))
    error = error/(param.nx * param.ny)
    return error
end

function solver(param)
    u = Array{Float64, 2}(undef, param.nx, param.ny)              # Initial solution
    unew = Array{Float64, 2}(undef, param.nx, param.ny)              # updated solution
    x = Array{Float64, 1}(undef, param.nx)             # grid in x direction
    y = Array{Float64, 1}(undef, param.ny)             # grid in y direction

    exact_sol = Array{Float64, 2}(undef, param.nx, param.ny)              # exact solution array
    
    # 2D grid generation Remove for loop use LinRange (check why is it efficient?)
    for i in 1:param.nx
        x = @. param.xmin + (i-1)*param.dx
    end
    exact_solution!(param, x, exact_sol)
    # Invoking initial condition
    initial_u!(param, x, u)

    t = 0.0
    it = 0.0

    dt = param.cfl/(abs(ux)/dx + abs(uy)/dy)
    sigma = abs(a) * dt / dx            # as a substitute to cfl


    while t < Tf 
        # -------Crucial block-------------
        if t + dt > Tf
            dt = Tf - t                    
            sigma = dt * abs(a) / dx
        end
        # ---------------------------------

        update_lw!(u, unew, sigma)
        u .= unew                       # Use . for element-wise operation on vector
        t += dt
        it += 1
    end
    err = error_cal(param, exact_sol, u)

    # Output to terminal    
    println("---------------------------")
    println("Error is: ", err)
    println("Iterations: ", it)
    println("---------------------------")

    # Writing solution to Files
    open("num_sol.txt","w") do io
        writedlm(io, [x u], "\t\t")
    end
    open("exact_sol.txt","w") do io
        writedlm(io, [x exact_sol], "\t\t")
    end
end

# for inputting parameters of the simulation
xmin, xmax = 0.0, 1.0                   # [xmin, xmax]
ymin, ymax = 0.0, 1.0                   # [ymin, ymax]

ux, uy = 1.0, 1.0                                  # velocity
nx, ny = 100, 100
Tf = 1                            # N = number of grid points, t = final time
cfl = 0.8
dx = (xmax - xmin)/(nx-1)
dy = (ymax - ymin)/(ny-1)

# @show N, Tf, xmin, xmax, 
param = (; N, Tf, dx, dy, xmin, ymin, ux, uy, cfl)
@timeit to "Solver" solver(param)
# show(to)