using Plots
using DelimitedFiles

num_data = readdlm("../linadv_ser/num_sol.txt", Float64)
exact_data = readdlm("../linadv_ser/exact_sol.txt", Float64)

plot(num_data[:,1],num_data[:,2], label="Exact Solution", linestyle=:solid, linewidth=2,dpi=150)
plot!(exact_data[:,1],exact_data[:,2],label="Numerical Solution", xlabel="Domain", ylabel="solution values(u)", title="Solution Plot",
    linewidth=2, linestyle=:dot, linecolor="black", dpi=150)

savefig("../linadv_ser/linadv_ser_plot.png")