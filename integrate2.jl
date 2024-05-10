#integration of cos(x) has been performed

using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

function integrate(a::Float64, b::Float64)
    return sin(b) - sin(a)
end

a, b = 0.0, 2.0
mya = a + ((b-a)/size) * rank
myb = mya + (b-a)/size

lres = integrate(mya, myb)
# if rank == 1
#     @show lres
# end
req = Vector{MPI.Request}(undef, size-1)
if rank == 0
    res = fill(0.0, size)
    res[1] = lres
    # @show res[1]
    req = Vector{MPI.Request}(undef, size-1)

    for i in 1:size-1
        req[i] = MPI.Irecv!(res, i, 0, comm)
    end
    @show res
    MPI.Waitall(req)
    @show res

    my_ans = 0.0
    for j in 1:size
        @show res
        global my_ans += res[j]
    end
    println("Result of integration is $my_ans")
else
    # if rank == 2
    #     @show lres
    # end
    MPI.send(lres, 0, 0, comm)
end

MPI.Finalize()


