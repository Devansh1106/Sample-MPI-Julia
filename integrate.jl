using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

# integrate for cos
function integrate(a::Float64, b::Float64)
    return b - a
end

leftB::Float64 = 0.0
rightB::Float64 = 2.0
step_size = (rightB - leftB)/size

left_ = *(step_size, rank)
right_ = +(step_size, left_)

res = integrate(left_, right_)
MPI.Barrier(comm)

recv_buf = MPI.Reduce(res, +, comm, root=0)

if rank == 0
    println("Result of the integration is: $(recv_buf)")
end