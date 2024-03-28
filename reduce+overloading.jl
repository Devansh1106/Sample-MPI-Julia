using MPI
MPI.Init()
import Base.+

struct Point{T}
    x::T
    y::T
end

+(A::Point{T}, B::Point{T}) where T = Point{T}(A.x + B.x, A.y + B.y)

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

p = Point(rank, rank)
print("Point at $(rank): $(p)\n")

recv_buf = MPI.Reduce(p, +, comm, root=0)
if rank == 0
    print("New point at rank 0: $(recv_buf)\n")
end
