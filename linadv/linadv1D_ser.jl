# solving linear advection equation:-
# u_t + a*u_x = 0 with periodic bc taking a = 1
# initial conditions: 
# u(x) = sin(2πx)
# Computational domain [0,1]

using DelimitedFiles
using Plots
using TimerOutputs
const to = TimerOutput()

# To generate a initial solution through initial condition
function initial_u!(param, x, u)
    for i in 1:param.N
        u[i] = sin(2.0 * π * x[i])
    end
end

# @views --> avoids unnecessary memory allocations while assignment by creating "view"
# a "view" provides a way to access and modify a subset of an array w/o creating a new array

# Lex-Wendroff method 
function update_lw!(u, unew, sigma)
    unew[1] = u[1] - 0.5* sigma *(u[2] - u[end-1]) + 0.5*sigma*sigma*(u[end-1] - 2*u[1] + u[2])   # u[0] = u[end-1] due to periodic bc
    @views unew[2:end-1] .= u[2:end-1] - 0.5*sigma*(u[3:end] - u[1:end-2]) + 0.5*sigma*sigma*(u[1:end-2] - 2*u[2:end-1] + u[3:end])
    unew[end] = u[end] - 0.5*sigma*(u[2] - u[end-1]) + 0.5*sigma*sigma*(u[end-1] - 2*u[end] + u[2]) # u[end+1] = u[2] due to periodic bc
end

# Exact solution calculation
function exact_solution!(param, x, exact_sol)
    for i in 1:param.N
        exact_sol[i] = sin(2.0 * π * (x[i]-(param.a * param.Tf)))
    end
end

# Error calculation
function error_cal(param, exact_sol, u)
    error = 0.0 
    error = sum(abs, exact_sol[1:param.N] - u[1:param.N])
    error = error/param.N
    return error
end


function solver(param)
    u = fill(0.0, param.N)              # Initial solution
    unew = fill(0.0, param.N)           # Updated solution
    x = fill(0.0, param.N)              # grid array
    exact_sol = fill(0.0, param.N)      # exact solution array
    
    # 1-D grid generation
    for i in 1:param.N
        x[i] = param.xmin + (i-1)*param.dx
    end
    exact_solution!(param, x, exact_sol)
    # Invoking initial condition
    initial_u!(param, x, u)

    t = 0.0
    it = 0.0

    dt = param.cfl * dx / abs(a)
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
a = 1                                   # velocity
N, Tf = 320, 1                          # N = number of grid points, t = final time
cfl = 0.8
dx = (xmax - xmin)/(N-1)
@show N, Tf, xmin, xmax, a, cfl
param = (; N, Tf, dx, xmin, a, cfl)
@timeit to "Solver" solver(param)
show(to)