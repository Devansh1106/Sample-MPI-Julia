using Plots
using DelimitedFiles

num_data = readdlm("numerical_parallel.txt", Float64)
exact_data = readdlm("exact_parallel.txt", Float64)

plot(num_data[:,1],num_data[:,2], 
     label="Exact Solution",
     linestyle=:solid, linewidth=2,
     dpi=150)

plot!(exact_data[:,1],exact_data[:,2], 
      label="Numerical Solution", xlabel="Domain", ylabel="solution values(u)",
      title="Solution Plot",
      linewidth=2, linestyle=:dot, linecolor="black", 
      dpi=150)

savefig("linadv1D_par.png")
println("Plot is in `linadv1D_par.png` file")