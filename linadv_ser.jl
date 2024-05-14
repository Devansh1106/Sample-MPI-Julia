# solving linear advection equation:-
# u_t + a*u_x = 0 with periodic bc taking a = 1
# initial conditions: 
# u(x) = sin(2πx)
# Computational domain [0,1]

# for inputting parameters of the simulation
function param_input()
    xmin, xmax = 0.0, 1.0  # [xmin, xmax]
    a = 1   # velocity
    N, t = 400, 1   # N= number of grid points, t = final time
    dx = (xmax - xmin)/(N-1)
    cfl = 0.8
    dt = (cfl * dx)/abs(a)
    @show dx
    @show dt
    @show N, t, cfl, xmin, xmax, a
    return N, t, dx, cfl, dt, xmin, a
end

# To generate a initial solution through initial condition
function initial_u!(N, x, u)
    for i in 1:N
        u[i] = sin(2.0 * π * x[i])
    end
    return u
end

# Lex-Wendroff method 
function update_lw!(u, cfl, unew)
    unew[1] = u[1] - 0.5*cfl*(u[2] - u[end-1]) + 0.5*cfl*cfl*(u[end-1] - 2*u[1] + u[2])   # u[0] = u[end-1] due to periodic bc
    unew[2:end-1] = u[2:end-1] - 0.5*cfl*(u[3:end] - u[1:end-2]) + 0.5*cfl*cfl*(u[1:end-2] - 2*u[2:end-1] + u[3:end])
    unew[end] = u[end] - 0.5*cfl*(u[2] - u[end-1]) + 0.5*cfl*cfl*(u[end-1] - 2*u[end] + u[2]) # u[end+1] = u[2] due to periodic bc
    return unew
end

# Exact solution calculation
function exact_solution!(x, N, t, a, exact_sol)
    for i in 1:N
        exact_sol[i] = sin(2.0 * π * (x[i]-(a * t)))
    end
    return exact_sol
end

# Error calculation
function error_cal(N, exact_sol, u)
    error = 0.0 
    error = sum(abs, exact_sol[1:N] - u[1:N])
    error = error/N
    return error
end


function solver()
    N, t, dx, cfl, dt, xmin, a = param_input()
    u = fill(0.0, N)        # Initial solution
    unew = fill(0.0, N)     # Updated solution
    x = fill(0.0, N)        # grid array
    exact_sol = fill(0.0, N) # exact solution array
    
    # 1-D grid generation
    for i in 1:N
        x[i] = xmin + (i-1)*dx
    end
    exact_sol_ = exact_solution!(x, N, t, a, exact_sol)

    # Invoking initial condition
    u_ini = initial_u!(N, x, u)

    j = 0.0
    it = 0.0

    while j <= t && it < 500
        # -------Crucial block-------------
        if j + dt > t
            dt = t - j                  # example: if j= 0.99 (<t) and dt = 0.5 hence if now loop runs it will give solution for final time
                                        # t = 0.99+0.5 which voilates the condition. Hence need to check this. Not very clear to me!!!
            cfl = dt * abs(a) / dx
        end
        # ---------------------------------

        unew = update_lw!(u_ini, cfl, unew)
        u_ini .= unew           # Use . for element-wise operation on vector
        j += dt
        it += 1
    end
    err = error_cal(N, exact_sol_, u)

    # Output to terminal
    println("---------------------------")
    println("Exact Solution Vector: ")
    println("---------------------------")
    println(exact_sol_[1])
    println(exact_sol_[2])
    println("⋮")
    println(exact_sol_[end-1])
    println(exact_sol_[end])
    
    println("---------------------------")
    println("Error is: ", err)
    println("Iterations: ", it)
    println("---------------------------")
    
    println("Numerical Solution vector: ")
    println("---------------------------")
    println(u[1])
    println(u[2])
    println(u[3])
    println("⋮")
    println(u[end-1])
    println(u[end])
end

solver()