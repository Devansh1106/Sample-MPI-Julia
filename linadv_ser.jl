# solving linear advection equation:-
# u_t + param.a*u_x = 0 with periodic bc taking param.a = 1
# initial conditions: 
# u(x) = sin(2πx)
# Computational domain [0,1]

# To generate param.a initial solution through initial condition
function initial_u!(param, x, u)
    for i in 1:param.N
        u[i] = sin(2.0 * π * x[i])
    end
end

# Lex-Wendroff method 
function update_lw!(u, unew, sigma)
    unew[1] = u[1] - 0.5* sigma *(u[2] - u[end-1]) + 0.5*sigma*sigma*(u[end-1] - 2*u[1] + u[2])   # u[0] = u[end-1] due to periodic bc
    @views unew[2:end-1] .= u[2:end-1] - 0.5*sigma*(u[3:end] - u[1:end-2]) + 0.5*sigma*sigma*(u[1:end-2] - 2*u[2:end-1] + u[3:end])
    unew[end] = u[end] - 0.5*sigma*(u[2] - u[end-1]) + 0.5*sigma*sigma*(u[end-1] - 2*u[end] + u[2]) # u[end+1] = u[2] due to periodic bc
end

# Exact solution calculation
function exact_solution!(x, param, exact_sol)
    for i in 1:param.N
        exact_sol[i] = sin(2.0 * π * (x[i]-(param.a * param.t)))
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
    u = fill(0.0, param.N)        # Initial solution
    unew = fill(0.0, param.N)     # Updated solution
    x = fill(0.0, param.N)        # grid array
    exact_sol = fill(0.0, param.N)   # exact solution array
    
    # 1-D grid generation
    for i in 1:param.N
        x[i] = param.xmin + (i-1)*param.dx
    end
    exact_solution!(x, param, exact_sol)

    # Invoking initial condition
    initial_u!(param, x, u)

    j = 0.0
    it = 0.0

    println("Enter the cfl: ")
    cfl = readline()
    cfl = parse(Float64, cfl)
    dt = cfl * dx / abs(a)
    sigma = abs(a) * dt / dx        # as a substitute to cfl


    while j <= t && it < 500
        # -------Crucial block-------------
        if j + dt > t
            dt = t - j                  # example: if j= 0.99 (<param.t) and param.dt = 0.5 hence if now loop runs it will give solution for final time
                                        # param.t = 0.99+0.5 which voilates the condition. Hence need to check this. Not very clear to me!!!
            sigma = dt * abs(a) / dx
        end
        # ---------------------------------

        update_lw!(u, unew, sigma)
        u .= unew           # Use . for element-wise operation on vector
        j += dt
        it += 1
    end
    err = error_cal(param, exact_sol, u)

    # Output to terminal    
    println("---------------------------")
    println("Error is: ", err)
    println("Iterations: ", it)
    println("---------------------------")
end

# for inputting parameters of the simulation
xmin, xmax = 0.0, 1.0  # [xmin, xmax]
a = 1   # velocity
N, t = 400, 1   # N= number of grid points, t = final time
dx = (xmax - xmin)/(N-1)
@show N, t, xmin, xmax, a
param = (; N, t, dx, xmin, a)
solver(param)