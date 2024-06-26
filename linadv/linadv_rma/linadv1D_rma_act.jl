# solving linear advection equation:-
# u_t + a*u_x = 0 with periodic bc taking a = 1
# initial conditions: 
# u(x) = sin(2πx)
# Computational domain [0,1]
# Exact solution: u(x-at) = sin(2π(x-at))

using DelimitedFiles
using Plots
using TimerOutputs
using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)
const to = TimerOutput()


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
Tf = parse(Int, ARGS[2])        # t = final time
cfl = parse(Float64, ARGS[3])
a = parse(Int, ARGS[4])         # velocity

xmin, xmax = 0.0, 1.0           # [xmin, xmax]
N_local = div(N, size)          # / -> converts N to float; div() doesn't
dx = (xmax - xmin)/(N - 1)

# if N % size !== 0 then last rank will have different number of cells
# otherwise, same.

# ----------------------------------------------------------------------
# Order of calculating `xmin_local` then `if...end` and then `xmax_local`
# should be maintained, it leads to correct distribution of grid when
# `N % size !== 0`
xmin_local = xmin + dx * rank * N_local
if rank == size - 1
    N_local = N - rank*N_local
end
xmax_local = xmin_local + dx * N_local
# -----------------------------------------------------------------------

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
function get_ghost_values!(param, u, win)
    if size == 1
        u[1] = u[param.N_local]
        u[param.N_local + 2] = u[3]
        return
    end
    next = (rank + 1) % size
    prev = (rank + size - 1) % size
    
    MPI.Win_fence(win)
    if rank !== size - 1  
        MPI.Put!(u[param.N_local+1], win; rank=next, disp=0)        # origin_buf, window; target_rank, disp
    else
        MPI.Put!(u[param.N_local], win; rank=next, disp=0)          # origin_buf, window, target_rank, disp
    end
    MPI.Win_fence(win)
    if rank == 0 
        disp_ = N - (size - 1)*N_local
        MPI.Put!(u[3], win; rank=prev, disp=disp_+1)                # origin_buf, window, target_rank, disp
    elseif rank !== 0 && rank !== size - 1
        MPI.Put!(u[2], win; rank=prev, disp=N_local + 1)            # origin_buf, window, target_rank, disp
    elseif rank == size -1
        disp_ = div(N, size)
        MPI.Put!(u[2], win; rank=prev, disp=disp_+1)                # origin_buf, window, target_rank, disp
    end
    MPI.Win_fence(win)
end

function collective_win_create(u)
    win = MPI.Win_create(u, comm)
    return win
end
    
function collective_win_free(win)
    MPI.free(win)
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
    err = fill(0.0, 1)

    # 1-D grid generation
    exact_solution!(param, x_local, exact_sol)

    # Invoking initial condition
    initial_u!(param, x_local, u)

    t = 0.0
    it = 0.0

    dt = param.cfl * param.dx / abs(param.a)
    sigma = abs(a) * dt / param.dx              # as a substitute to cfl

    win = collective_win_create(u)
    while t < Tf
        # -------Crucial block-------------
        if t + dt > Tf
            dt = Tf - t                         # if `j + dt` goes beyond `t` then loop will quit in next iteration hence solution will not get calculated till `t`
                                                # With this, solution will get calculated as close to `t`
            sigma = dt * abs(param.a) / param.dx
        end
        # ---------------------------------
        get_ghost_values!(param, u, win)

        update_lw!(u, unew, sigma)
        @views u[2:end-1] .= unew[2:end-1]      # Use . for element-wise operation on vector
        t += dt
        it += 1
    end
    collective_win_free(win)
    win_err = collective_win_create(err)
    err[1] = error_cal(param, exact_sol, u)
    MPI.Win_fence(win_err)
    MPI.Accumulate!(err[1], MPI.SUM, win_err; rank=0, disp=0)
    collective_win_free(win_err)
    
    if rank == 0
        # Output to terminal    
        println("---------------------------")
        println("Error is: ", err)
        println("Iterations: ", it)
        println("---------------------------")
    end

    # Writing solution to Files
    open("num_sol_par_$rank.txt", "w") do io
        writedlm(io, [x_local u[2:end-1]], "\t\t")
    end

    open("exact_sol_par_$rank.txt", "w") do io
        writedlm(io, [x_local exact_sol], "\t\t")
    end
end
@timeit to "Solver" solver(param)
if rank == 0
    show(to)
end