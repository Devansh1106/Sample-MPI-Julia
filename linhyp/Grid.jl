module Grid

struct CartesianGrid
    domain::Tuple{Float64, Float64} # xmin, xmax
    nx::Int64                   # nx - number of points
    xc::Vector{Float64}         # cell centers
    xf::Vector{Float64}         # cell faces
    dx::Vector{Float64}         # cell sizes
end

# Uniform cartesian grid 
function make_grid(problem::Problem, param::Parameters)
    nvar = problem.nvar
    xmin, xmax = problem.domain 
    nx = param.grid_size
    println("Making Uniform Grid of interval [$xmin, $xmax]")
    dx = (xmax - xmin)/nx
    xc = LinRange(xmin+0.5*dx, xmax-0.5*dx, nx)
    println("Grid with number of points = $nx")
    println("(xmin, xmax) = ($xmin, $xmax)")
    xf = LinRange(xmin, xmax, nx+1)
    return CartesianGrid((xmin, xmax), nx, xc, xf, dx)
end

export make_grid

end # module