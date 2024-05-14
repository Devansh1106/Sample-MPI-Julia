# solving linear advection equation:-
# u_t + a*u_x = 0 with periodic bc taking a = 1
# initial conditions: 
# u(x) = sin(2πx)
# Computational domain [0,1]

function solver()
    xmin, xmax = 0.0, 1.0  # [xmin, xmax]
    a = 1   # velocity
    N, t = 5, 1   # N= number of grid points, t = final time
    dx = (xmax - xmin)/(N-1)
    cfl = 0.8
    dt = (cfl * dx)/abs(a)
    # println(dx)
    return N, t, dx, cfl, dt, xmin, a
end

function initial_u(N, x, u)
    for i in 1:N
        u[i] = sin(2.0 * π * x[i])
    end
    # println(u[1])
    # println(u[2])
    # println(u[end-1])
    # println(u[end])
    # @show u
    return u
end

function update_lw(u, cfl, unew)
    unew[1] = u[1] - 0.5*cfl*(u[2] - u[end-1]) + 0.5*cfl*cfl*(u[end-1] - 2*u[1] + u[2])   # u[0] = u[end-1] due to periodic bc
    unew[2:end-1] = u[2:end-1] - 0.5*cfl*(u[3:end] - u[1:end-2]) + 0.5*cfl*cfl*(u[1:end-2] - 2*u[2:end-1] + u[3:end])
    unew[end] = u[end] - 0.5*cfl*(u[2] - u[end-1]) + 0.5*cfl*cfl*(u[end-1] - 2*u[end] + u[2]) # u[end+1] = u[2] due to periodic bc

    return unew
end

N, t, dx, cfl, dt, xmin, a= solver()
u = fill(0.0, N)
unew = fill(0.0, N)
x = fill(0.0, N)
# t_ = fill(0.0, N)
for i in 1:N
    global x[i] = xmin + (i-1)*dx
    # global t_[i] = (i-1)*dt
end
@show x
u = initial_u(N, x, u)
j = 0
it = 0

while j < t && it < 150
    global unew = update_lw(u, cfl, unew)
    global u = unew
    global j += dt
    global it += 1
end

function exact_solution(x, N, t, a)
    exact_sol = fill(0.0, N)
    for i in 1:N
        exact_sol[i] = sin(2.0 * π * (x[i]-(a * t)))
    end
    return exact_sol
end

exact_sol = exact_solution(x, N, t, a)
println("Exact Solution Vector: ")
# println(exact_sol[1])
# println(exact_sol[2])
# println("⋮")
# println(exact_sol[end-1])
# println(exact_sol[end])
@show exact_sol
error = 0.0 
error = sum(abs, exact_sol[1:N] - u[1:N])
error = error/N

println("Error is: ",error)
println("Iterations: ", it)
println("Numerical Solution vector: ")
@show u
# println(u[1])
# println(u[2])
# println(u[3])
# println("⋮")
# println(u[end-1])
# println(u[end])

