module HypSys1D

using LinearAlgebra
using OffsetArrays
using DelimitedFiles
using Plots
using Grid
using EqLinAdv

struct Problem{F <: Function}
    domain::Tuple{Float64, Float64}
    nvar::Int64
    initial_value::F
    # boundary_value::F
    boundary_condition::Any
    final_time::Float64
end

struct Scheme
    equation
    numflux::String
end

struct Parameters
    cfl::Float64
    grid_size::Int64
end

function create_parameters(cfl::Float64, grid_size::Int64)
    @assert (cfl >= 0.0) "cfl should be >= 0.0"
    return Parameters(cfl, grid_size)
end

gArray(nvar::Int64, nx::Int64) = OffsetArray(zeros(nvar, nx+2),
                                             OffsetArrays.Origin(1,0))

function adjust_time_step(problem::Problem, param::Parameters,
                          dt::Float64, t::Float64)
    final_time = problem.final_time
    if t + dt > final_time
        dt = final_time - t
    end

    return dt    
end

# TODO: for non-uniform grid there will be loop and dt = min(...);lam=max(...)
function compute_lam_dt(equation, param::Parameters, scheme::Scheme, grid::CartesianGrid)
    if scheme.numflux == "rusanov"
        lam = maximum(abs.(eigvals(equation.fprime)))
    end
    dt = (param.cfl * grid.dx[1])/lam   # dx is a vector; it is useful for non-uniform grid version
    return lam, dt
end

function set_initial_value!(grid::CartesianGrid, U::gArray, problem::Problem) where gArray
    nx = grid.nx
    xc = grid.xc
    initial_value = problem.initial_value
    for i in 1:nx
        @views U[:,i] .= initial_value(xc[i])
    end
    # println(U)
end

function update_ghost!(grid::CartesianGrid, U::gArray, problem::Problem) where gArray
    nx = grid.nx
    if problem.boundary_condition == "Periodic"
        U[:, 0] .= U[:, nx-1]
        U[:, nx+1] .= U[:, 2]
    end
    # display(U)
end

function compute_residual!(equation, grid::CartesianGrid, lam::Float64, U::gArray,
                           scheme::Scheme, res::gArray, problem::Problem) where gArray
    Uf = zeros(problem.nvar)
    nx = grid.nx
    # xf = grid.xf
    dx = OffsetArray(zeros(nx+2), OffsetArrays.Origin(0))
    @. dx[1:nx] = grid.dx[1:nx]
    dx[0] = grid.dx[nx-1]
    dx[nx+1] = grid.dx[2]
    numflux = scheme.numflux
    for i in 1:nx+1
        @views Ul, Ur = U[:, i-1], U[:, i]
        # display(Ur)
        if scheme.numflux == "rusanov"
            rusanov!(equation, lam, Ul, Ur, Uf)
        end
        # display(Uf)
        @views res[:, i-1] += Uf/ dx[i-1]           # How is this working for non-uniform grid??
        @views res[:, i] -= Uf/ dx[i]
    end
    # display(res)
end

function solve(equation, problem::Problem, scheme::Scheme, param::Parameters)
    grid = make_grid(problem, param)
    nvar = problem.nvar
    Tf = problem.final_time
    nx = grid.nx
    dx = grid.dx
    xf = grid.xf
    # Allocating variables
    U = gArray(nvar, nx+1)
    Ue = gArray(nvar, nx)
    res = gArray(nvar, nx) # dU/dt + res(U) = 0
    set_initial_value!(grid, U, problem)
    # display(U)
    it, t = 0, 0.0
    while t < Tf
        lam, dt = compute_lam_dt(equation, param, scheme, grid)
        # @show (lam, dt)
        dt = adjust_time_step(problem, param, dt, t)
        # @show (lam, dt)
        update_ghost!(grid, U, problem)
        # display(U)
        compute_residual!(equation, grid, lam, U, scheme, res, problem)
        @. U -= dt*res
        # display(U)
        t += dt; it += 1
        # update_plot!(grid, problem, equation, scheme, U, t, it, param, plt_data)
    end
    # @show U
    compute_exact_soln!(grid, equation, problem, t, Ue)
    # println(U,"\n")
    # println(Ue)
    # plot_sol(grid, U, Ue, problem)

    # open("linhyp/num_sol.txt", "w") do io
    #     writedlm(io, [U], '\n')
    # end
    # println("Solution is in num_sol.txt file")
end

function plot_sol(grid::CartesianGrid, U::gArray, Ue::gArray, problem::Problem) where gArray
    nvar = problem.nvar
    nx = grid.nx
    xc = grid.xc
    xf = grid.xf
    for i in 1:nvar
        @views y_exact, y_num = Ue[i, 1:nx], U[i, 1:nx+1]
        plot(xc, y_exact, label="Exact Solution", linestyle=:solid, linewidth=2, dpi=150)
        plot!(xf, y_num, label="Numerical Solution", xlabel="Domain", ylabel="Solution values",
              title="Solution Plot", linewidth=2, linestyle=:dot, linecolor="black", 
              dpi=150)
        savefig("linhyp/fig/hypsys1D_$i.png")
    end
    println("Figures are saved in folder `linhyp/fig`.")
end

export Problem
export Parameters
export create_parameters
export Scheme
export solve
export gArray

end # module