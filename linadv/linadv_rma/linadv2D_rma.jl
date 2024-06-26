# solving linear advection equation 2D:-
# u_t + a*u_x + b*u_y= 0 with periodic bc taking a = b = 1
# initial conditions: 
# u(xy) = sin(2πx) * sin(2πy)
# Computational domain [0,1] × [0,1]

using DelimitedFiles
using TimerOutputs
using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
nprocs = MPI.Comm_size(comm)
const to = TimerOutput()

# To generate a initial solution through initial condition
function initial_u!(param, x, y, u)
    for i in 2:param.nx+1, j in 2:param.ny+1
        u[i,j] = sin(2.0 * π * x[i-1]) * sin(2.0 * π * y[j-1])
    end
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
function error_cal!(param, exact_sol, u)
    err = 0.0
    for i in 2:param.nx+1, j in 2:param.ny+1
        err += abs(exact_sol[i-1,j-1] - u[i,j])
    end
    err = err/(param.nx * param._ny)
    return err
end

function halo_exchange!(param, u, win)
    # Updating top and bottom row
    for j in 2:param.ny+1
        u[1,j] = u[end-1,j] 
        u[end,j] = u[2,j]
    end

    if nprocs == 1        
        # Updating left and right columns
        for i in 1:param.nx+2
            u[i,1] = u[i,end-1] 
            u[i,end] = u[i,2]
        end
        return
    else
        next = (rank + 1) % nprocs
        prev = (rank + nprocs - 1) % nprocs
        buf = fill(0.0, param.nx+2)
        buf1 = fill(0.0, param.nx+2)

        # Sending the second last column of each rank to `next` rank
        for i in 1:param.nx+2
            buf[i] = u[i, end-1]
        end
        MPI.Win_lock(win; rank=next, type=MPI.LOCK_SHARED, nocheck=true)    # window; target_rank, lock_type, nocheck::bool
        MPI.Put!(buf, win; rank=next, disp=0)                               # origin_buf, window; target_rank, disp
        MPI.Win_unlock(win, rank=next)

        # Sending the second column of each rank to `prev` rank
        for i in 1:param.nx+2
            buf1[i] = u[i,2]
        end

        # Displacement calculation according to target_rank
        if rank == 0
            d = param._ny - param.ny * (nprocs-1)                           # param._ny is undistributed grid size in y 
            _disp = (d+1) * (param.nx+2)
        elseif rank == nprocs - 1
            d = div(param._ny, nprocs)
            _disp = (d+1) * (param.nx+2)
        else
            _disp = (param.ny+1) * (param.nx+2)
        end
        MPI.Win_lock(win; rank=prev, type=MPI.LOCK_SHARED, nocheck=true)    # window; target_rank, lock_type, nocheck::bool
        MPI.Put!(buf1, win; rank=prev, disp=_disp)                          # origin_buf, window; target_rank, disp
        MPI.Win_unlock(win, rank=prev)
    end
    MPI.Win_fence(win)
end

# Window creation routine
function collective_win_create(u)
    win = MPI.Win_create(u, comm)
    return win
end
    
# Window free routine
function collective_win_free(win)
    MPI.free(win)
end

function solver(param)
    u = Array{Float64, 2}(undef, param.nx+2, param.ny+2)       # Initial solution
    unew = Array{Float64, 2}(undef, param.nx, param.ny)        # updated solution
    x = Array{Float64, 1}(undef, param.nx)                     # grid in x direction
    y = Array{Float64, 1}(undef, param.ny)                     # grid in y direction
    exact_sol = Array{Float64, 2}(undef, param.nx, param.ny)   # exact solution array
    err = fill(0.0, 1)
    
    # 2D grid generation
    x = LinRange(param.xmin+0.5*param.dx, param.xmax-0.5*param.dx, nx)
    y = LinRange(param.ymin+0.5*param.dy, param.ymax-0.5*param.dy, ny)

    # Exact solution
    @timeit to "exact_solution!" begin
        exact_sol = exact_solution!(param, x, y, exact_sol)
    end

    # Invoking initial condition
    @timeit to "initial_u!" begin
        initial_u!(param, x, y, u)
    end

    t = 0.0
    it = 0.0

    dt = param.cfl/(abs(a)/dx + abs(b)/dy)
    sigma_x = abs(a) * dt / dx            # as a substitute to cfl
    sigma_y = abs(b) * dt / dy            # as a substitute to cfl

    win = collective_win_create(u)
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
            halo_exchange!(param, u, win)
        end

        @timeit to "update_lw!" begin
            update_lw!(param, u, unew, sigma_x, sigma_y)
        end

        for i in 1:param.nx, j in 1:param.ny
            u[i+1, j+1] = unew[i, j]        # Use . for element-wise operation on vector
        end
        t += dt
        it += 1
    end
    collective_win_free(win)                # free the window of `u` array

    # Window creation for error calculation
    win_err = collective_win_create(err)

    @timeit to "error_cal!" begin
        err[1] = error_cal!(param, exact_sol, u)
    end
    MPI.Win_fence(win_err)
    MPI.Accumulate!(err[1], MPI.SUM, win_err; rank=0, disp=0)

    collective_win_free(win_err)            # free the error window 

    if rank == 0
        # Output to terminal    
        println("---------------------------")
        println("Error is: ", err[1])
        println("Iterations: ", it)
        println("---------------------------")
        
        # Writing undistributed part of mesh (here x) in file
        open("mesh_x.txt", "w") do io
            writedlm(io, x)
        end
    end
    # Writing solution to files
    open("num_sol2D_par_$rank.txt","w") do io
        writedlm(io, u[2:end-1,2:end-1])
    end
    open("exact_sol2D_par_$rank.txt","w") do io
        writedlm(io, exact_sol)
    end
    open("mesh_y_$rank.txt","w") do io
        writedlm(io, y)
    end
end

@timeit to "Global variables" begin
    # for inputting parameters of the simulation
    xmin, xmax = 0.0, 1.0                 # [xmin, xmax]
    _ymin, _ymax = 0.0, 1.0               # [ymin, ymax]

    a, b = 1.0, 1.0                       # velocity
    nx, _ny = 200, 200
    Tf = 1
    cfl = 0.8
    dx = (xmax - xmin)/(nx)
    dy = (_ymax - _ymin)/(_ny)

    ny = div(_ny, nprocs)                 # Dividing `ny` among processes

    # if `_ny` is not divisible by `nprocs` then 
    # `ny` will have `remainder` columns 
    if rank == nprocs - 1
        ny = _ny - (ny * rank)
    end

    # Distributing domain in `y` among `nprocs`
    ymin = _ymin + (rank*dy*ny)
    ymax = ymin + (dy*ny)

    param = (; nx, ny, _ny, Tf, dx, dy, xmin, xmax, ymin, ymax, a, b, cfl)
    @show param
end

@timeit to "Solver" solver(param)
if rank == 0
    show(to)
end