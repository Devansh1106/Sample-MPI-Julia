module EqLinAdv

using LinearAlgebra
using DelimitedFiles
using Grid
using HypSys1D

struct LinAdv
    fprime::Array{Float64, 2}
end

function flux(fprime::Array{Float64, 2}, U::SubArray{Float64, 1})
    F = fprime * U
    # display(U)
    return F
end

function rusanov!(equation::LinAdv, lam::Float64, Ul::SubArray{Float64, 1},
                  Ur::SubArray{Float64, 1}, Uf::Vector{Float64})
    fprime = equation.fprime
    Fl, Fr = flux(fprime, Ul), flux(fprime, Ur)
    Uf .= 0.5*(Fl + Fr) - 0.5*lam*(Ur - Ul)
    # display(Uf)
end

# TODO: Check expression for exact solution from notes.
function compute_exact_soln!(grid::CartesianGrid, equation::LinAdv, problem::Problem,
                             t::Float64, Ue::gArray) where gArray
    nx = grid.nx
    xc = grid.xc
    nvar = problem.nvar
    initial_value = problem.initial_value
    eigen_vals = eigvals(equation.fprime)
    eigen_vecs = eigvecs(equation.fprime)
    for j in 1:nx
        for i in 1:nvar
            Ue[i,j] = (inv(eigen_vecs) * initial_value(xc[j] - eigen_vals[i] * t))[i]
        end
    end
    for i in 1:nx
        @views Ue[:,i] .= eigen_vecs * Ue[:, i]
    end
    # display(Ue)
    # open("linhyp/exact_sol.txt", "w") do io
    #     writedlm(io, [Ue], '\n')
    # end
    # println("Exact Solution is in exact_sol.txt")
end

export rusanov!
export compute_exact_soln!
export LinAdv

end # module