using Plots
using DelimitedFiles
plotlyjs()
num_data = readdlm("num_sol2D.txt", Float64)  
exact_data = readdlm("exact_sol2D.txt", Float64)

plot(num_data[:,1],num_data[:,2],num_data[:,3:end],
     label="Numerical Solution", st=:surface, xlabel="x", ylabel="y",
     title="Solution Plot", zlabel="Numerical solution",dpi=150)

savefig("linadv2D_ser_num.html")

plot(exact_data[:,1],exact_data[:,2], exact_data[:,3:end],
      label="Exact Solution", xlabel="x", ylabel="y", zlabel="Exact solution",
      title="Solution Plot", st=:surface, dpi=150)

savefig("linadv2D_ser_exact.html")
println("Plot is in `linadv2D_ser_exact.html` and `linadv2D_ser_num.html` file")