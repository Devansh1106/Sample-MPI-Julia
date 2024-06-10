using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

buf = fill(0, 1)
buf[1] = rank

win = MPI.Win_create(buf, comm)
MPI.Win_fence(win)
MPI.Accumulate!(buf, MPI.SUM, win; rank=0, disp=0)
MPI.Win_fence(win)
@show buf
if rank == 0
    println(buf)
end