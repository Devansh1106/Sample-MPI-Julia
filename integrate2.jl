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
if rank == 0
    my_ans = 0.0
    res = fill(0.0, size)
    res[1] = lres
    # @show res
    req = Vector{MPI.Request}(undef, size-1)
    for i in 1:size-1
        MPI.Irecv!(res[i+1], i, comm, req[i])
        # println("ji")
    end
    MPI.Waitall!(req)
    for num in res
        global my_ans += num
    end
    print("Result of integration is $my_ans")
else
    MPI.send(lres, comm, dest=0)
end

MPI.Finalize()


