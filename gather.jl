using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

send_buf = nothing
if rank == 0
    send_buf = collect(1:size) .* 100
    print("Original array on rank 0: $(send_buf)\n")
end

recv_buf = MPI.Scatter(send_buf, Int, comm, root=0)
recv_buf *= recv_buf
v = MPI.Gather(recv_buf, comm, root=0)
if rank == 0
    print("New array after using gather operation:\n $(v)\n")
end
