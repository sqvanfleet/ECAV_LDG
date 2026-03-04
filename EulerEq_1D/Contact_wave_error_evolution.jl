using StartUpDG
using LinearAlgebra
using JLD2
using StaticArrays
using Plots
using LaTeXStrings

base_path = @__DIR__

data = load(joinpath(base_path, "Contact_wave_data", "N2_K80_tspan4.0.jld2"), "rd", "md", "sol", "L2_error")

rd = data[1]
md = data[2]
sol = data[3]
L2_error = data[4]


common_size = (800, 600)
margins = (5Plots.mm, 5Plots.mm, 5Plots.mm, 5Plots.mm) 

p1 = plot(sol.t[2:end], L2_error[2:end],
        yaxis = :log10, 
        size = common_size,
        left_margin = margins[1],
        right_margin = margins[2],
        top_margin = margins[3],
        bottom_margin = margins[4],
        #ylims = (-14,2),
        xlabel = L"$t$", 
        #ylabel = L"$\log_{10}||u_h - u||_{L^2(D)}$", 
        label = "N = $(rd.N)" , linewidth = 3,
        tickfontsize = 20,
        legendfontsize = 16,
        xguidefontsize = 24,)


for i in [4,6,8]
    local loop_file = joinpath(base_path, "Contact_wave_data", "N$(i)_K80_tspan4.0.jld2")
    local data = load(loop_file, "rd", "md", "sol", "L2_error")
    local rd = data[1]
    local md = data[2]
    local sol = data[3]
    local L2_error = data[4]

    plot!(p1,sol.t[2:end], L2_error[2:end], 
        label = "N = $(rd.N)", linewidth = 3)

end

savefig(p1, "Contact_wave_plots/L2_error_evolution.pdf")
