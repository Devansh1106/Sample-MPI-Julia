using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

send_buf = nothing
recv_buf = Vector{Float64}(undef, size)

if rank == 0
    send_buf = rand(Float64, (size,size))
    print("Original array on rank 0:\n $(send_buf)\n")
end

MPI.Scatter!(send_buf, recv_buf, comm, root=0) # recv_buf will be given inside MPI.Scatter!
print("I got this array on $(rank):\n$(recv_buf)\n")
