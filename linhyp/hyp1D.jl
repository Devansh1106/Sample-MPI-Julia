module HypSys1D

using OffsetArrays
using Grid

struct Problem{F <: Function}
    domain::Tuple{Float64, Float64}
    nvar::Int64
    initial_value::F
    boundary_value::F
    boundary_condition::Any
    final_time::Float64
end

struct Parameters
    cfl::Float64
    grid_size::Int64
end

function Parameters(cfl::Float64, grid_size::Int64)
    @assert (cfl >= 0.0) "cfl should be >= 0.0"
    return Parameters(cfl, grid_size)
end

struct Scheme
    equation::String
    numflux::String
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

# TODO: Need clarification for lambda calculation
function compute_lam_dt(equation, scheme::Scheme, param::Parameters,
                        grid::Grid, )

end

function set_initial_value!(grid::CartesianGrid, U::Vector{Float64}, problem::Problem)
    nx = grid.nx
    xc = grid.xc
    initial_value = problem.initial_value
    for i in nx
        U[:,i] = initial_value(xc[i])
    end
end

function update_ghost!(grid::CartesianGrid, U::Vector{Float64}, problem::Problem)
    xmin, xmax = grid.domain
    nx = grid.nx
    if problem.boundary_condition == "Periodic"
        U[:, 0] .= U[:, nx-1]
        U[:, nx+1] .= U[:, 2]
    end
end

function compute_residual!(grid::CartesianGrid, lam::Float64, U::Vector{Float64},
                           scheme::Scheme, res::Float64)
    Uf = zeros(grid.nvar)
    nx = grid.nx
    xf = grid.xf
    dx = grid.dx
    numflux = scheme.numflux
    for i in 1:nx+1
        
end



end # module