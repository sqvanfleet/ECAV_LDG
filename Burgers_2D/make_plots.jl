using LaTeXStrings
using Plots
using JLD2
using OrdinaryDiffEq
using StartUpDG

common_size = (800, 600)
margins = (5Plots.mm, 5Plots.mm, 5Plots.mm, 5Plots.mm) # (left, right, top, bottom)

data = load("Data/N3_K30_LDG_modal_EC.jld2", "md", "rd", "sol")

sol_LDG = data[3]

data = load("Data/N3_K30_BR1_modal_EC.jld2", "md", "rd", "sol")

sol_BR1 = data[3]

dSdt_evolution_plot = plot(sol_LDG.prob.p.t, sol_LDG.prob.p.dSdt,
                size = common_size,
                left_margin = margins[1],
                right_margin = margins[2],
                top_margin = margins[3],
                bottom_margin = margins[4],
                xlabel = L"$t$", 
                xticks = [0.0,0.1,0.3,0.5,0.7],
                #ylabel = L"$\frac{\mathrm{d}}{\mathrm{d}t}\int_{D}S$",
                #title = "modal",
                #yguidefontrotation = 0, 
                label = "LDG",
                linewidth = 3, tickfontsize = 20, legendfontsize = 16,
                xguidefontsize = 24)
plot!(sol_BR1.prob.p.t, sol_BR1.prob.p.dSdt,
                label = "BR1", linewidth = 3)


# 1. Create the "Plots" directory (does nothing if it already exists)
mkpath("Plots")

# 2. Save the plot
savefig(dSdt_evolution_plot,"Plots/dSdt_plot_Burgers_2D.png")

epsilon_evolution_plot = plot(sol_LDG.prob.p.t, sol_LDG.prob.p.max_epsilon .+ 1e-14,
                yscale = :log10,
                size = common_size,
                left_margin = margins[1],
                right_margin = margins[2],
                top_margin = margins[3],
                bottom_margin = margins[4],
                xlabel = L"$t$", 
                xticks = [0.0,0.1,0.3,0.5,0.7],
                #ylabel = L"\epsilon",
                #title = "modal",
                #yguidefontrotation = 0, 
                label = "LDG", linewidth = 3, tickfontsize = 20, legendfontsize = 16,
                xguidefontsize = 24,
                legend =:topleft)
plot!(sol_BR1.prob.p.t, sol_BR1.prob.p.max_epsilon .+ 1e-14,
                label = "BR1", linewidth = 3)

savefig(epsilon_evolution_plot,"Plots/LDG_epsilon_evolution.png")

print("LDG used $(sol_LDG.destats.naccept) time steps, BR1 used $(sol_BR1.destats.naccept) time steps")
# display(epsilon_evolution_plot)

# data = load("Data/N3_K30_LDG_nodal_EC.jld2", "md", "rd", "sol")

# sol_LDG = data[3]

# data = load("Data/N3_K30_BR1_nodal_EC.jld2", "md", "rd", "sol")

# sol_BR1 = data[3]

# dSdt_evolution_plot = plot(sol_LDG.prob.p.t, sol_LDG.prob.p.dSdt,
#                 size = common_size,
#                 left_margin = margins[1],
#                 right_margin = margins[2],
#                 top_margin = margins[3],
#                 bottom_margin = margins[4],
#                 xlabel = L"$t$", 
#                 ylabel = L"$\frac{\mathrm{d}}{\mathrm{d}t}\int_{D}S$",
#                 title = "nodal",
#                 yguidefontrotation = 0, 
#                 label = "LDG time steps = $(sol_LDG.destats.naccept)", linewidth = 3, tickfontsize = 20, legendfontsize = 16,
#                 xguidefontsize = 24)
# plot!(sol_BR1.prob.p.t, sol_BR1.prob.p.dSdt,
#                 label = "BR1 time steps = $(sol_BR1.destats.naccept)", linewidth = 3)

# display(dSdt_evolution_plot)

# epsilon_evolution_plot = plot(sol_LDG.prob.p.t, sol_LDG.prob.p.max_epsilon,
#                 size = common_size,
#                 left_margin = margins[1],
#                 right_margin = margins[2],
#                 top_margin = margins[3],
#                 bottom_margin = margins[4],
#                 xlabel = L"$t$", 
#                 ylabel = L"\epsilon",
#                 title = "nodal",
#                 yguidefontrotation = 0, 
#                 label = "LDGtime steps = $(sol_LDG.destats.naccept)", linewidth = 3, tickfontsize = 20, legendfontsize = 16,
#                 xguidefontsize = 24)
# plot!(sol_BR1.prob.p.t, sol_BR1.prob.p.max_epsilon,
#                 label = "BR1time steps = $(sol_BR1.destats.naccept)", linewidth = 3)


# display(epsilon_evolution_plot)