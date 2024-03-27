using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)
println(rank)

# println("Hello world $(rank) out of $(size)")
buf::Vector{Int64} = zeros(1)
if rank == 0
    buf = [1]
    MPI.Send(buf, 1, 10, comm)
end
if rank == 1
    MPI.Recv!(buf, 0, 10, comm)  #MPI.Recv(buf, 0, 10, comm)---> Will not work. MPI.Recv is for receiving a single `isbits` objects
    println(buf)
end