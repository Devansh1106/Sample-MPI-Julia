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
    for i in 2:param.N_local + 1
        u[i] = sin(2.0 * π * x[i-1])
    end
end

# getting boundary values and halo exchanges
function get_ghost_values!(param, u)
    if size == 1
        u[1] = u[param.N_local]
        u[param.N_local + 2] = u[3]
        return
    end
    next = (rank + 1) % size
    prev = (rank + size - 1) % size
    
    # Things matters while communicating:
    # 1- Order of sending messages
    # 2- tag matching
    
    if rank == 0
        MPI.send(u[N_local+1], comm, dest=next, tag=0)
        MPI.send(u[3], comm, dest=prev, tag=0)
        u[N_local + 2] = MPI.recv(comm; source=next, tag=0)
        u[1] = MPI.recv(comm; source=prev, tag=0)
    elseif rank == size - 1
        MPI.send(u[2], comm, dest=prev, tag=0)
        MPI.send(u[N_local], comm, dest=next, tag=0)
        u[1] = MPI.recv(comm; source=prev, tag=0)
        u[N_local + 2] = MPI.recv(comm; source=next, tag=0)
    else
        MPI.send(u[N_local + 1], comm, dest=next, tag=0)
        MPI.send(u[2], comm, dest=prev, tag=0)
        u[N_local + 2] = MPI.recv(comm; source=next, tag=0)
        u[1] = MPI.recv(comm; source=prev, tag=0)
    end
    # @show (rank, u)
    # exit()
end


# Lex-Wendroff method 
function update_lw!(u, unew, sigma)
    @views unew[2:end-1] .= u[2:end-1] - 0.5*sigma*(u[3:end] - u[1:end-2]) + 0.5*sigma*sigma*(u[1:end-2] - 2.0*u[2:end-1] + u[3:end])
end

# Exact solution calculation
function exact_solution!(param, x, exact_sol)
    for i in 1:param.N_local
        exact_sol[i] = sin(2.0 * π * (x[i]-(param.a * param.t)))
    end
end

# Error calculation
function error_cal(param, exact_sol, u)
    error = 0.0 
    error = sum(abs, exact_sol[1:param.N_local] - u[2:param.N_local+1])
    error = error/param.N_local/size
    return error
end

# for inputting parameters of the simulation
xmin, xmax = 0.0, 1.0           # [xmin, xmax]
a = 1                           # velocity
N, t = 100, 1                   # N = number of grid points, t = final time
N_local = div(N, size)          # / -> converts N to float; div() doesn't;   change in error_cal also
x = nothing
dx = (xmax - xmin)/(N - 1)
if rank == 0
    x = fill(0.0, N)            # grid array
    for i in 1:N
        x[i] = xmin + (i-1)*dx
    end
    # @show x
    # exit()
end

x_local = fill(0.0, N_local)
MPI.Scatter!(x, x_local, comm, root=0)

grid_per_rank = (xmax - xmin)/size
xmin = rank*grid_per_rank
xmax = xmin + grid_per_rank
# if rank == 1
#     @show xmin, xmax
#     exit()
# end

@show N_local, t, xmin, xmax, a, dx
param = (; N_local, t, dx, xmin, a)


function solver(param)
    u = fill(0.0, param.N_local + 2)            # Initial solution
    unew = fill(0.0, param.N_local + 2)         # Updated solution
    exact_sol = fill(0.0, param.N_local)        # exact solution array
    
    # 1-D grid generation
    exact_solution!(param, x_local, exact_sol)

    # Invoking initial condition
    initial_u!(param, x_local, u)

    j = 0.0
    it = 0.0
    buf = fill(0.0, 2)

    if rank == 0
        print("Enter the cfl: ")
        cfl = readline()                # should be less than 1.0
        cfl = parse(Float64, cfl)       # parse() will convert it to Float64
        dt = cfl * dx / abs(a)          # readline() input it as string
        sigma = abs(a) * dt / dx        # as a substitute to cflu[2]
        buf[1] = dt
        buf[2] = sigma
    end
    buf = MPI.bcast(buf, 0, comm)


    while j <= t && it < 500
        # -------Crucial block-------------
        if j + buf[1] > t
            buf[1] = t - j                  # example: if j= 0.99 (<param.t) and param.dt = 0.5 hence if now loop runs it will give solution for final time
                                            # param.t = 0.99+0.5 which voilates the condition. Hence need to check this. Not very clear to me!!!
            buf[2] = buf[1] * abs(a) / dx
        end
        # ---------------------------------

        get_ghost_values!(param, u)
        # @show it
        # if it == 0
        #     @show (rank, u)
        #     exit()
        # end

        update_lw!(u, unew, buf[2])
        @views u[2:end-1] .= unew[2:end-1]      # Use . for element-wise operation on vector
        j += buf[1]
        it += 1
    end
    local_err = error_cal(param, exact_sol, u)
    global_err = MPI.Reduce(local_err, +, comm, root=0)

    if rank == 0
        final_sol = fill(0.0, size * N_local)
        x_final = fill(0.0, size * N_local)
        exact_sol_final = fill(0.0, size * N_local)
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

    # u_final = fill(0.0, param.N_local)
    u_final = @view u[2:end-1]
    # if rank == 0
    #     @show u
    #     exit()
    # end
    final_sol = MPI.Gather(u_final, comm, root=0)
    x_final = MPI.Gather(x_local, comm, root=0)
    exact_sol_final = MPI.Gather(exact_sol, comm, root=0)


    
    # Writing solution to Files
    open("/home/devansh/Sample-MPI-Julia/linadv_ser/num_sol_par.txt","w") do io
        writedlm(io, [x_final final_sol], "\t\t")
    end
    open("/home/devansh/Sample-MPI-Julia/linadv_ser/exact_sol_par.txt","w") do io
        writedlm(io, [x_final exact_sol_final], "\t\t")
    end

    # # Plotting: saved in "linadv_ser.png"
    # plot(x_local, exact_sol, label="Exact Solution", linestyle=:solid, linewidth=2,dpi=150)
    # plot!(x_local, u, label="Numerical Solution", xlabel="Domain", ylabel="solution values(u)", title="Solution Plot",
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