using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

send_buf = nothing

if rank == 0
    send_buf = collect(1:size) .* 100
    print("Original array on rank 0:\n $(send_buf)\n")
end

recv_buf = MPI.Scatter(send_buf, Int, comm, root=0)  # recv_buf = MPI.Scatter(buffer(all data), datatype, ::COMM, root::Int64(rank with all data))
print("I got this on rank $(rank):\n $(recv_buf)\n")