module EqLinAdv

struct LinAdv{nvar <: Int64}
    fprime::Array{nvar, nvar}
end
# function get_equation(fprime::Array{nvar, nvar}) where nvar <: Int64

# end


end # module