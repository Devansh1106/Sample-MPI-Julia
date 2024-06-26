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
    for i in 2:param.nx+1, j in 2:param.ny+1
        u[i,j] = sin(2.0 * π * x[i-1]) * sin(2.0 * π * y[j-1])
    end
    return u
end

# Lex-Wendroff method 
function update_lw!(param, u, unew, sigma_x, sigma_y)
    for i in 2:param.nx+1, j in 2:param.ny+1
        unew[i-1,j-1] = (u[i,j] - 0.5 * sigma_x * (u[i+1,j] - u[i-1,j])
                                - 0.5 * sigma_y * (u[i,j+1] - u[i,j-1])
                                + 0.5 * sigma_x^2 * (u[i-1,j] - 2.0*u[i,j] + u[i+1,j])
                                + 0.25 * sigma_x * sigma_y * (u[i+1,j+1] - u[i-1,j+1] - u[i+1,j-1] + u[i-1,j-1])
                                + 0.5 * sigma_y^2 * (u[i,j-1] - 2.0*u[i,j] + u[i,j+1]))
    end
end

# Exact solution calculation
function exact_solution!(param, x, y, exact_sol)
    _a = param.a
    _b = param.b
    _Tf = param.Tf
    for i in 1:param.nx, j in 1:param.ny
        exact_sol[i,j] = sin(2.0 * π * (x[i]-(_a * _Tf))) * sin(2.0 * π * (y[j]-(_b * _Tf)))
    end
    return exact_sol
end

# Error calculation
function error_cal!(param, err, exact_sol, u)
    for i in 2:param.nx+1, j in 2:param.ny+1
        err += abs(exact_sol[i-1,j-1] - u[i,j])
    end
    err = err/(param.nx * param.ny)
    return err
end

function halo_exchange!(param, u)
    # Updating top and bottom row
    for j in 2:param.ny+1
        u[1,j] = u[end-1,j] 
        u[end,j] = u[2,j]
    end
    
    # Updating left and right columns
    for i in 1:param.nx+2
        u[i,1] = u[i,end-1] 
        u[i,end] = u[i,2]
    end
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
    @timeit to "exact_solution!" begin
        exact_sol = exact_solution!(param, x, y, exact_sol)
    end

    # Invoking initial condition
    @timeit to "initial_u!" begin
        u = initial_u!(param, x, y, u)
    end

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
        @timeit to "halo_exchange!" begin
            halo_exchange!(param, u)
        end

        @timeit to "update_lw!" begin
            update_lw!(param, u, unew, sigma_x, sigma_y)
        end
        @views u[2:end-1,2:end-1] .= unew        # Use . for element-wise operation on vector
        t += dt
        it += 1
    end
    @timeit to "error_cal!" begin
        err = error_cal!(param, err, exact_sol, u)
    end

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
dx = (xmax - xmin)/nx
dy = (ymax - ymin)/ny

param = (; nx, ny, Tf, dx, dy, xmin, xmax, ymin, ymax, a, b, cfl)
@show param
@timeit to "Solver" solver(param)
show(to)