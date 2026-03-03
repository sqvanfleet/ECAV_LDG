using JLD2

N_1_K_8 = load("Density_Wave_convergence_data/N1_K8_tspan10.0_nodal_BR1.jld2", "rd", "md", "sol")
N_1_K_16 = load("Density_Wave_convergence_data/N1_K16_tspan10.0_nodal_BR1.jld2", "rd", "md", "sol")
N_1_K_32 = load("Density_Wave_convergence_data/N1_K32_tspan10.0_nodal_BR1.jld2", "rd", "md", "sol")
N_1_K_64 = load("Density_Wave_convergence_data/N1_K64_tspan10.0_nodal_BR1.jld2", "rd", "md", "sol")

sol_8 = N_1_K_8[3]
sol_16 = N_1_K_16[3]
sol_32 = N_1_K_32[3]
sol_64 = N_1_K_64[3]

sol_8_error = sol_8.prob.p.L2_error[end]
sol_16_error = sol_16.prob.p.L2_error[end]
sol_32_error = sol_32.prob.p.L2_error[end]
sol_64_error = sol_64.prob.p.L2_error[end]

errors_N_1 = [sol_8_error,
              sol_16_error,
              sol_32_error,
              sol_64_error]

error_ratio_N_1 = (1/log(2))*log.([
    sol_8_error / sol_16_error,
    sol_16_error / sol_32_error,
    sol_32_error / sol_64_error
])

N_2_K_8 = load("Density_Wave_convergence_data/N2_K8_tspan10.0_nodal_BR1.jld2", "rd", "md", "sol")
N_2_K_16 = load("Density_Wave_convergence_data/N2_K16_tspan10.0_nodal_BR1.jld2", "rd", "md", "sol")
N_2_K_32 = load("Density_Wave_convergence_data/N2_K32_tspan10.0_nodal_BR1.jld2", "rd", "md", "sol")
N_2_K_64 = load("Density_Wave_convergence_data/N2_K64_tspan10.0_nodal_BR1.jld2", "rd", "md", "sol")

sol_8 = N_2_K_8[3]
sol_16 = N_2_K_16[3]
sol_32 = N_2_K_32[3]
sol_64 = N_2_K_64[3]

sol_8_error = sol_8.prob.p.L2_error[end]
sol_16_error = sol_16.prob.p.L2_error[end]
sol_32_error = sol_32.prob.p.L2_error[end]
sol_64_error = sol_64.prob.p.L2_error[end]

errors_N_2 = [sol_8_error,
              sol_16_error,
              sol_32_error,
              sol_64_error]

error_ratio_N_2 = (1/log(2))*log.([
    sol_8_error / sol_16_error,
    sol_16_error / sol_32_error,
    sol_32_error / sol_64_error
])

N_3_K_8 = load("Density_Wave_convergence_data/N3_K8_tspan10.0_nodal_BR1.jld2", "rd", "md", "sol")
N_3_K_16 = load("Density_Wave_convergence_data/N3_K16_tspan10.0_nodal_BR1.jld2", "rd", "md", "sol")
N_3_K_32 = load("Density_Wave_convergence_data/N3_K32_tspan10.0_nodal_BR1.jld2", "rd", "md", "sol")
N_3_K_64 = load("Density_Wave_convergence_data/N3_K64_tspan10.0_nodal_BR1.jld2", "rd", "md", "sol")

sol_8 = N_3_K_8[3]
sol_16 = N_3_K_16[3]
sol_32 = N_3_K_32[3]
sol_64 = N_3_K_64[3]

sol_8_error = sol_8.prob.p.L2_error[end]
sol_16_error = sol_16.prob.p.L2_error[end]
sol_32_error = sol_32.prob.p.L2_error[end]
sol_64_error = sol_64.prob.p.L2_error[end]

errors_N_3 = [sol_8_error,
              sol_16_error,
              sol_32_error,
              sol_64_error]

error_ratio_N_3 = (1/log(2))*log.([
    sol_8_error / sol_16_error,
    sol_16_error / sol_32_error,
    sol_32_error / sol_64_error
])


N_4_K_8 = load("Density_Wave_convergence_data/N4_K8_tspan10.0_nodal_BR1.jld2", "rd", "md", "sol")
N_4_K_16 = load("Density_Wave_convergence_data/N4_K16_tspan10.0_nodal_BR1.jld2", "rd", "md", "sol")
N_4_K_32 = load("Density_Wave_convergence_data/N4_K32_tspan10.0_nodal_BR1.jld2", "rd", "md", "sol")
N_4_K_64 = load("Density_Wave_convergence_data/N4_K64_tspan10.0_nodal_BR1.jld2", "rd", "md", "sol")

sol_8 = N_4_K_8[3]
sol_16 = N_4_K_16[3]
sol_32 = N_4_K_32[3]
sol_64 = N_4_K_64[3]

sol_8_error = sol_8.prob.p.L2_error[end]
sol_16_error = sol_16.prob.p.L2_error[end]
sol_32_error = sol_32.prob.p.L2_error[end]
sol_64_error = sol_64.prob.p.L2_error[end]

errors_N_4 = [sol_8_error,
              sol_16_error,
              sol_32_error,
              sol_64_error]

error_ratio_N_4 = (1/log(2))*log.([
    sol_8_error / sol_16_error,
    sol_16_error / sol_32_error,
    sol_32_error / sol_64_error
])

# display(error_ratio_N_1)
# display(error_ratio_N_2)
# display(error_ratio_N_3)
# display(error_ratio_N_4)

# display(errors_N_1)
# display(errors_N_2)
# display(errors_N_3)
# display(errors_N_4)

using Plots

sol_32_epsilon = sol_32.prob.p.max_epsilon
sol_32_t = sol_32.prob.p.t

plot(sol_32_t, sol_32_epsilon, xlabel="Time", ylabel="Max Epsilon", title="Max Epsilon vs Time for N=3, K=32", legend=false)


using Printf

K_values = [8, 16, 32, 64]

# Store all errors and ratios in dictionaries for easy looping
errors = Dict(
    1 => errors_N_1,
    2 => errors_N_2,
    3 => errors_N_3,
    4 => errors_N_4
)

ratios = Dict(
    1 => error_ratio_N_1,
    2 => error_ratio_N_2,
    3 => error_ratio_N_3,
    4 => error_ratio_N_4
)

# Initialize table string
table_str = ""
table_str *= "\\begin{table}[h]\n"
table_str *= "\\centering\n"
table_str *= "\\begin{tabular}{|c | c c | c c | c c | c c|}\n"
table_str *= "\\hline\n"
table_str *= "K & \\multicolumn{2}{|c|}{N=1} & \\multicolumn{2}{|c|}{N=2} & \\multicolumn{2}{|c|}{N=3} & \\multicolumn{2}{c|}{N=4} \\\\\n"
table_str *= "  & Error & Rate & Error & Rate & Error & Rate & Error & Rate \\\\\n"
table_str *= "\\hline\n"

# Fill rows
for (i, K) in enumerate(K_values)
    row = @sprintf("%d", K)
    for N in 1:4
        err = errors[N][i]
        # Rate should be "-" for the first row (no previous data)
        rate = i == 1 ? NaN : ratios[N][i - 1]
        row *= @sprintf(" & %.3e & %s", err, isnan(rate) ? "-" : @sprintf("%.2f", rate))
    end
    global table_str *= row * " \\\\\n"
end

table_str *= "\\hline\n\\end{tabular}\n\\caption{Convergence study errors and rates}\n\\end{table}\n"

println(table_str)

