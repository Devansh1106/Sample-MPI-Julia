# solving linear advection equation:-
# u_t + a*u_x = 0 with periodic bc taking a = 1
# initial conditions: 
# u(x) = sin(2πx)
# Computational domain [0,1]

using DelimitedFiles
using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)


# To generate a initial solution through initial condition
function initial_u!(param, x, u)
    for i in 1:param.N
        u[i] = sin(2.0 * π * x[i])
    end
end

# Lex-Wendroff method 
function update_lw!(u, unew, sigma)
    # local_start = rank*N + 1
    # local_end = (rank + 1) * N
    # unew[local_start] = u[local_start] - 0.5* sigma *(u[local_start + 1] - u[local_end-1]) + 0.5*sigma*sigma*(u[local_end-1] - 2*u[local_start] + u[local_start+1])   # u[0] = u[end-1] due to periodic bc
    # @views unew[local_start + 1:local_end - 1] .= u[local_start+1:local_end-1] - 0.5*sigma*(u[local_start+2:local_end] - u[local_start:local_end-2]) + 0.5*sigma*sigma*(u[local_start:local_end-2] - 2*u[local_start+1:local_end-1] + u[local_start+2:local_end])
    # unew[local_end] = u[local_end] - 0.5*sigma*(u[local_start+1] - u[local_end-1]) + 0.5*sigma*sigma*(u[local_end-1] - 2*u[local_end] + u[local_start+1]) # u[end+1] = u[2] due to periodic bc
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
    error = error/param.N/size
    return error
end

# for inputting parameters of the simulation
xmin, xmax = 0.0, 1.0  # [xmin, xmax]
a = 1   # velocity
N, t = 100, 1   # N = number of grid points, t = final time
N = div(N, size)    # / -> converts N to float; div() doesn't
grid_per_rank = (xmax - xmin)/size
xmin = 0.0 + rank*grid_per_rank
xmax = xmin + grid_per_rank
dx = (xmax - xmin)/(N - grid_per_rank)

@show N, t, xmin, xmax, a, dx
param = (; N, t, dx, xmin, a)


function solver(param)
    u = fill(0.0, param.N)        # Initial solution
    unew = fill(0.0, param.N)     # Updated solution
    x = fill(0.0, param.N)        # grid array
    exact_sol = fill(0.0, param.N)   # exact solution array
    
    # 1-D grid generation
    for i in 1:param.N
        x[i] = param.xmin + (i-1)*param.dx
    end
    exact_solution!(param, x, exact_sol)
    # Invoking initial condition
    initial_u!(param, x, u)

    j = 0.0
    it = 0.0
    buf = fill(0.0, 2)

    if rank == 0
        println("Enter the cfl: ")
        cfl = readline()                # should be less than 1.0
        cfl = parse(Float64, cfl)       # parse() will convert it to Float64
        dt = cfl * dx / abs(a)          # readline() input it as string
        sigma = abs(a) * dt / dx        # as a substitute to cfl
        # buf = (dt, sigma)
        buf[1] = dt
        buf[2] = sigma
        # MPI.bcast(buf, comm, root=0)
    end
    # println(buf)
    # buf::Tuple{Float64, Float64}
    buf = MPI.bcast(buf, 0, comm)



    while j <= t && it < 500
        # -------Crucial block-------------
        if j + buf[1] > t
            buf[1] = t - j                  # example: if j= 0.99 (<param.t) and param.dt = 0.5 hence if now loop runs it will give solution for final time
                                        # param.t = 0.99+0.5 which voilates the condition. Hence need to check this. Not very clear to me!!!
            buf[2] = buf[1] * abs(a) / dx
        end
        # ---------------------------------

        update_lw!(u, unew, buf[2])
        u .= unew           # Use . for element-wise operation on vector
        j += buf[1]
        it += 1
    end
    local_err = error_cal(param, exact_sol, u)
    global_err = MPI.Reduce(local_err, +, comm, root=0)

    if rank == 0
        final_sol = fill(0.0, size * N)
        x_final = fill(0.0, size * N)
        exact_sol_final = fill(0.0, size * N)
        # Output to terminal    
        println("---------------------------")
        println("Error is: ", global_err)
        println("Iterations: ", it)
        println("---------------------------")
    end
    # # for inputting parameters of the simulation
    # xmin, xmax = 0.0, 1.0  # [xmin, xmax]
    # a = 1   # velocity
    # N, t = 100, 1   # N = number of grid points, t = final time
    # dx = (xmax - xmin)/(N-1)
    # @show N, t, xmin, xmax, a
    # param = (; N, t, dx, xmin, a)


    final_sol = MPI.Gather(u, comm, root=0)
    x_final = MPI.Gather(x, comm, root=0)
    exact_sol_final = MPI.Gather(exact_sol, comm, root=0)


    
    # Writing solution to Files
    open("/home/devansh/Sample-MPI-Julia/linadv_ser/num_sol_par.txt","w") do io
        writedlm(io, [x_final final_sol], "\t\t")
    end
    open("/home/devansh/Sample-MPI-Julia/linadv_ser/exact_sol_par.txt","w") do io
        writedlm(io, [x_final exact_sol_final], "\t\t")
    end

    # # Plotting: saved in "linadv_ser.png"
    # plot(x, exact_sol, label="Exact Solution", linestyle=:solid, linewidth=2,dpi=150)
    # plot!(x, u, label="Numerical Solution", xlabel="Domain", ylabel="solution values(u)", title="Solution Plot",
    #     linewidth=2, linestyle=:dot, linecolor="black", dpi=150)
    # savefig("../linadv_ser/linadv_ser.png")
end

# # for inputting parameters of the simulation
# xmin, xmax = 0.0, 1.0  # [xmin, xmax]
# a = 1   # velocity
# N, t = 100, 1   # N = number of grid points, t = final time
# dx = (xmax - xmin)/(N-1)
# @show N, t, xmin, xmax, a
# param = (; N, t, dx, xmin, a)
solver(param)
# solver(param)