# solving linear advection equation:-
# u_t + a*u_x = 0 with periodic bc taking a = 1
# initial conditions: 
# u(x) = sin(2πx)
# Computational domain [0,1]

using DelimitedFiles

# To generate a initial solution through initial condition
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
function exact_solution!(param, x, exact_sol)
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
