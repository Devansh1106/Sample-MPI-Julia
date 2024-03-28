using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

if rank == 0
    buf = [1,2]
    print("Original is: $(buf)\n")
    MPI.send(buf, comm, dest=1)
end
if rank == 1
    buf = MPI.recv(comm, source=0)  #MPI.Recv(buf, 0, 10, comm)---> Will not work. MPI.Recv is for receiving a single `isbits` objects
    println(buf)
end