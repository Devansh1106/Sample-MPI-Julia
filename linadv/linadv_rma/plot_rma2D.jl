using Plots
using DelimitedFiles
plotlyjs()

num_data = readdlm(joinpath(@__DIR__,"numerical_parallel.txt"), Float64)  # @__DIR__ - this macro expnands to a string with directory path of the file
exact_data = readdlm(joinpath(@__DIR__,"exact_parallel.txt"), Float64)

plot(num_data[:,1],num_data[:,2],num_data[:,3:end],
     label="Numerical Solution", st=:surface, xlabel="x", ylabel="y",
     title="Solution Plot", zlabel="Numerical solution",dpi=150)
savefig("linadv2D_par_num.html")

plot(exact_data[:,1],exact_data[:,2], exact_data[:,3:end],
      label="Exact Solution", xlabel="x", ylabel="y", zlabel="Exact solution",
      title="Solution Plot", st=:surface, dpi=150)

savefig("linadv2D_par_exact.html")
println("Plot is in `linadv2D_par_exact.html` and `linadv2D_par_num.html` file")