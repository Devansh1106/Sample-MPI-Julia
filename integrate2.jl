# integration of cos(x) has been performed

# Not fully correct. Working on the issue...

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
    @show req                                  
    for i in 1:size-1
        req[i] = MPI.Irecv!(res, i, 0, comm)
    end
    @show req
    @show res
else
    MPI.send(lres, 0, 0, comm)
end

if rank == 0
    my_ans = lres   # *
    # validate_requests(req)
    # MPI.isnull(req)
    # @show req
    MPI.Waitall(req)    # Waiting should be before the line in which communicated data need to be used
    for j in 1:size-1
        global my_ans += res[j] # global otherwise loop's scope will create it own my_ans variable instead of using * one
    end
    println("Result of integration is $my_ans")
end

# MPI.Finalize()    # Not needed in julia (also in tmux it will create a problem in running the program more than once in a single session)