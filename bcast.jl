using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

buf = nothing
if rank == 0
    buf = collect(1:5)
    print("This is the bcast array:\n $(buf)\n")
end

recv_buf = MPI.bcast(buf, comm, root=0)
print("I am rank $(rank) and got $(recv_buf)\n")
