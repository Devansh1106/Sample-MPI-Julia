using MPI
MPI.Init()

comm = MPI.COMM_WORLD
println("Hello World, rank$(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm))")
MPI.Barrier(comm)