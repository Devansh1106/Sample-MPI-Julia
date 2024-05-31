using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

function integrate(a::Float64, b::Float64)
    return sin(b) - sin(a)
    # return b - a
end

a, b = 0.0, 1.0
mya = a + ((b-a)/size) * rank
myb = mya + (b-a)/size

lres = integrate(mya, myb)

res = MPI.Gather(lres, comm, root=0)

if rank == 0
    my_ans = 0.0
    for j in 1:size
        global my_ans += res[j]
    end
    println("Result of integration is $my_ans")
end
