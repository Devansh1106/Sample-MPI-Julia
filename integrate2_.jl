using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

function integrate(a::Float64, b::Float64)
    return sin(b) - sin(a)
end

a, b = 0.0, 15.0
mya = a + ((b-a)/size) * rank
myb = mya + (b-a)/size

lres = integrate(mya, myb)

if rank == 0
    res = fill(0.0, size-1)
    req = Vector{MPI.Request}(undef, size-1)

    for i in 1:size-1
        req[i] = MPI.Irecv!(res, i, 0, comm)
    end
    
    MPI.Waitall(req) # Wait for all data to be received before continuing

    my_ans = lres
    for j in 1:size-1
        global my_ans += res[j]
    end
    println("Result of integration is $my_ans")
else
    MPI.send(lres, 0, 0, comm)
end

# MPI.Finalize()
