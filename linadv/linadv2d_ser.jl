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
    @views u[2:end-1,2:end-1] = @. sin(2.0 * π * x) * sin(2.0 * π * y)
    return u
end

# Lex-Wendroff method 
function update_lw!(u, unew, sigma_x, sigma_y)
    @views unew .= (u[2:end-1,2:end-1] - 0.5 * sigma_x * (u[3:end,2:end-1] - u[1:end-2,2:end-1])
                                       - 0.5 * sigma_y * (u[2:end-1,3:end] - u[2:end-1,1:end-2])
                                       + 0.5 * sigma_x^2 * (u[1:end-2,2:end-1] - 2.0*u[2:end-1,2:end-1] + u[3:end,2:end-1])
                                       + 0.25 * sigma_x * sigma_y * (u[3:end,3:end] - u[1:end-2,3:end] - u[3:end,1:end-2] + u[1:end-2,1:end-2])
                                       + 0.5 * sigma_y^2 * (u[2:end-1,1:end-2] - 2.0*u[2:end-1,2:end-1] + u[2:end-1,3:end]))
    return unew
end

# Exact solution calculation
function exact_solution!(param, x, y, exact_sol)
    _a = param.a
    _b = param.b
    _Tf = param.Tf
    exact_sol = @. sin(2.0 * π * (x-(_a * _Tf))) * sin(2.0 * π * (y-(_b * _Tf)))
    return exact_sol
end

# Error calculation
function error_cal!(param, err, exact_sol, u)
    @views err = sum(abs.(exact_sol[1:end,1:end] - u[2:end-1,2:end-1]))
    err = err/(param.nx * param.ny)
    return err
end

function halo_exchange!(u, x, y)
    # Updating top and bottom row
    @views u[1,2:end-1] .= u[end-1,2:end-1] 
    @views u[end,2:end-1] .= u[2,2:end-1]
    
    # Updating left and right columns
    @views u[2:end-1,1] .= u[2:end-1,end-1] 
    @views u[2:end-1,end] .= u[2:end-1,2]

    # Updating corners
    u[1,1] = u[end-1,end-1]
    u[end,end] = u[2,2]
    u[end,1] = u[2,end-1]
    u[1,end] = u[end-1, 2]
    return u
end

function solver(param)
    u = Array{Float64, 2}(undef, param.nx+2, param.ny+2)       # Initial solution
    unew = Array{Float64, 2}(undef, param.nx, param.ny)        # updated solution
    x = Array{Float64, 1}(undef, param.nx)                     # grid in x direction
    y = Array{Float64, 1}(undef, param.ny)                     # grid in y direction
    exact_sol = Array{Float64, 2}(undef, param.nx, param.ny)   # exact solution array
    err = 0.0
    
    # 2D grid generation
    x = LinRange(param.xmin+0.5*param.dx, param.xmax-0.5*param.dx, nx)
    y = LinRange(param.ymin+0.5*param.dy, param.ymax-0.5*param.dy, ny)
    x = x'

    # Exact solution
    exact_sol = exact_solution!(param, x, y, exact_sol)

    # Invoking initial condition
    u = initial_u!(param, x, y, u)

    t = 0.0
    it = 0.0

    dt = param.cfl/(abs(a)/dx + abs(b)/dy)
    sigma_x = abs(a) * dt / dx            # as a substitute to cfl
    sigma_y = abs(b) * dt / dy            # as a substitute to cfl


    while t < Tf 
        # -------Crucial block-------------
        if t + dt > Tf
            dt = Tf - t                    
            sigma_x = dt * abs(param.a) / param.dx
            sigma_y = dt * abs(param.b) / param.dy
        end
        # ---------------------------------

        # halo exchange
        u = halo_exchange!(u, x, y)
        unew = update_lw!(u, unew, sigma_x, sigma_y)
        @views u[2:end-1,2:end-1] .= unew        # Use . for element-wise operation on vector
        t += dt
        it += 1
    end
    err = error_cal!(param, err, exact_sol, u)

    # Output to terminal    
    println("---------------------------")
    println("Error is: ", err)
    println("Iterations: ", it)
    println("---------------------------")

    # Writing solution to Files     How are we going to write this in files?
    open("num_sol2D.txt","w") do io
        writedlm(io, [x' y u[2:end-1,2:end-1]], "\t\t")
    end
    open("exact_sol2D.txt","w") do io
        writedlm(io, [x' y exact_sol], "\t\t")
    end
end

# for inputting parameters of the simulation
xmin, xmax = 0.0, 1.0                   # [xmin, xmax]
ymin, ymax = 0.0, 1.0                   # [ymin, ymax]

a, b = 1.0, 1.0                         # velocity
nx, ny = 200, 200
Tf = 1
cfl = 0.8
dx = (xmax - xmin)/(nx)
dy = (ymax - ymin)/(ny)

param = (; nx, ny, Tf, dx, dy, xmin, xmax, ymin, ymax, a, b, cfl)
@show param
@timeit to "Solver" solver(param)
show(to)