# solving linear advection equation:-
# u_t + a*u_x = 0 with periodic bc taking a = 1
# initial conditions: 
# u(x) = sin(2πx)
# Computational domain [0,1]
# Exact solution: u(x-at) = sin(2π(x-at))

using DelimitedFiles
using Plots
using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)
rank == 0 ? time_start = MPI.Wtime() : nothing

# for inputting parameters of the simulation
if rank == 0
    if length(ARGS) != 4
        println("ArgumentError: Usage: julia linadv1D_par.jl <N> <t> <cfl> <a>")
        exit()
    end
else
    if length(ARGS) != 4
        exit()
    end
end

N = parse(Int, ARGS[1])         # N = number of grid points
@assert N > size
Tf = parse(Int, ARGS[2])         # t = final time
cfl = parse(Float64, ARGS[3])
a = parse(Int, ARGS[4])         # velocity

xmin, xmax = 0.0, 1.0           # [xmin, xmax]
N_local = div(N, size)          # / -> converts N to float; div() doesn't
dx = (xmax - xmin)/(N - 1)

# if N % size !== 0 then last rank will have different number of cells
# otherwise, same.
xmin_local = xmin + dx * rank * N_local
if rank == size - 1
    N_local = N - rank*N_local
end
xmax_local = xmin_local + dx * N_local

x_local = fill(0.0, N_local)
for i in 1:N_local
    x_local[i] = xmin_local + (i-1)*dx
end

@show N_local, Tf, xmin_local, xmax_local, a, cfl, dx
param = (; N_local, Tf, dx, a, cfl)


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
end


# Lex-Wendroff method 
function update_lw!(u, unew, sigma)
    @views unew[2:end-1] .= u[2:end-1] - 0.5*sigma*(u[3:end] - u[1:end-2]) + 0.5*sigma*sigma*(u[1:end-2] - 2.0*u[2:end-1] + u[3:end])
end

# Exact solution calculation
function exact_solution!(param, x, exact_sol)
    for i in 1:param.N_local
        exact_sol[i] = sin(2.0 * π * (x[i]-(param.a * param.Tf)))
    end
end

# Error calculation
function error_cal(param, exact_sol, u)
    error = 0.0 
    error = sum(abs, exact_sol[1:param.N_local] - u[2:param.N_local+1])
    error = error/param.N_local/size
    return error
end

function solver(param)
    u = fill(0.0, param.N_local + 2)            # Initial solution
    unew = fill(0.0, param.N_local + 2)         # Updated solution
    exact_sol = fill(0.0, param.N_local)        # exact solution array
    
    # 1-D grid generation
    exact_solution!(param, x_local, exact_sol)

    # Invoking initial condition
    initial_u!(param, x_local, u)

    t = 0.0
    it = 0.0

    dt = param.cfl * param.dx / abs(param.a)
    sigma = abs(a) * dt / param.dx        # as a substitute to cfl

    while t < Tf
        # -------Crucial block-------------
        if t + dt > Tf
            dt = Tf - t                      # if `j + dt` goes beyond `t` then loop will quit in next iteration hence solution will not get calculated till `t`
                                            # With this, solution will get calculated as close to `t`
            sigma = dt * abs(param.a) / param.dx
        end
        # ---------------------------------

        get_ghost_values!(param, u)

        update_lw!(u, unew, sigma)
        @views u[2:end-1] .= unew[2:end-1]      # Use . for element-wise operation on vector
        t += dt
        it += 1
    end
    local_err = error_cal(param, exact_sol, u)
    global_err = MPI.Reduce(local_err, +, comm, root=0)

    if rank == 0
        # Output to terminal    
        println("---------------------------")
        println("Error is: ", global_err)
        println("Iterations: ", it)
        println("---------------------------")
    end

    rank == 0 ? time_end = MPI.Wtime() : nothing
    rank == 0 ? println("Time taken: $(time_end - time_start)") : nothing

    # Writing solution to Files
    open("../linadv/num_sol_par_$rank.txt", "w") do io
        writedlm(io, [x_local u[2:end-1]], "\t\t")
    end

    open("../linadv/exact_sol_par_$rank.txt", "w") do io
        writedlm(io, [x_local exact_sol], "\t\t")
    end

    if rank == 0
        # Plotting: saved as "linadv1D_par.png"
        run(`sh -c "cat num_sol_par_*.txt > numerical_parallel.txt"`)
        run(`sh -c "cat exact_sol_par_*.txt > exact_parallel.txt"`)

        run(`sh -c "rm num_sol_par_*.txt"`)
        run(`sh -c "rm exact_sol_par_*.txt"`)

        num_data = readdlm("numerical_parallel.txt", Float64)
        exact_data = readdlm("exact_parallel.txt", Float64)

        plot(num_data[:,1],num_data[:,2], 
             label="Exact Solution",
             linestyle=:solid, linewidth=2,
             dpi=150)

        plot!(exact_data[:,1],exact_data[:,2], 
              label="Numerical Solution", xlabel="Domain", ylabel="solution values(u)",
              title="Solution Plot",
              linewidth=2, linestyle=:dot, linecolor="black", 
              dpi=150)

        savefig("linadv1D_par.png")
        println("DONE")
    end
end

solver(param)