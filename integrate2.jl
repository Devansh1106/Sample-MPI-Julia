#integration of cos(x) has been performed

# Currently due to a bug, code will give correct answer if used with n=1 process only.
# Running on n>1 will return error.
# Working on the issue...

using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

function integrate(a::Float64, b::Float64)
    return sin(b) - sin(a)
end

a, b = 0.0, 16.0
mya = a + ((b-a)/size) * rank # # Check for errors
# for i in 1:size-1
#     if MPI.Test_cancelled(statuses[i])
#         println("Request $i was cancelled")
#     elseif MPI.Test_status(statuses[i])
#         error_code = MPI.Error_code(statuses[i])
#         println("Request $i completed with error code $error_code")
#     else
#         println("Request $i completed successfully")
#     end
# end
myb = mya + (b-a)/size

lres = integrate(mya, myb)

if rank == 0
    res = fill(0.0, size-1)
    req = Vector{MPI.Request}(undef, size-1)

    for i in 1:size-1
        req[i] = MPI.Irecv!(res, i, 0, comm)
    end
else
    MPI.send(lres, 0, 0, comm)
end

if rank == 0
    my_ans = lres   # *
    @show res
    statuses = MPI.Waitall!(req)
    @show statuses    # Waiting should be before the line in which communicated data need to be used

    # # Check for errors
    # for i in 1:size-1
    #     if MPI.Test_cancelled(statuses[i])
    #         println("Request $i was cancelled")
    #     elseif MPI.Test_status(statuses[i])
    #         error_code = MPI.Error_code(statuses[i])
    #         println("Request $i completed with error code $error_code")
    #     else
    #         println("Request $i completed successfully")
    #     end
    # end

    for j in 1:size-1
        global my_ans += res[j] # global otherwise loop's scope will create it own my_ans variable instead of using * one
    end
    println("Result of integration is $my_ans")
end

# MPI.Finalize()    # Not needed in julia (also in tmux it will create a problem in running the program more than once in a single session)