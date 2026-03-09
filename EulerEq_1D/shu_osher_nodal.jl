using StartUpDG
using LinearAlgebra
using OrdinaryDiffEq, RecursiveArrayTools
using JLD2
using Trixi
using Trixi.ForwardDiff
using StaticArrays
using Plots
using LaTeXStrings

initial_condition_type = :shu_osher

include("initial_conditions.jl")
include("rhs.jl")

#Choose AV type
# AV_type = :BR1
AV_type = :LDG

#Choose modal or nodal DG
#DG_type = :modal
DG_type = :nodal

gamma = 1.4
equations = CompressibleEulerEquations1D(gamma)

N = 3
K1D = 100 
if DG_type == :modal
    rd = RefElemData(Line(), N) #modal
elseif DG_type == :nodal
    rd = RefElemData(Line(),SBP(), N) #nodal
end

(VX,), EToV = uniform_mesh(rd.element_type, K1D)

@. VX = 5.0*VX
md = MeshData((VX,), EToV, rd)
u = rd.Pq * shu_osher.(md.xq, equations)

tspan = (0.0, 1.8)

folder = joinpath(@__DIR__, "Shu_Osher_data")
mkpath(folder)
AV_type = :LDG
t = Float64[]
max_epsilon = Float64[]
du_visc_norm = Float64[]
sigma_norm = Float64[]
du_Euler_norm = Float64[]
du_Euler_AV_norm = Float64[]
params = (; rd, md, equations, interface_flux = flux_hllc, 
    initial_condition_type, AV_type, t, max_epsilon, DG_type,
    du_visc_norm, sigma_norm, du_Euler_norm, du_Euler_AV_norm)
filename = "N$(N)_K$(K1D)_tspan$(tspan[2])_$(DG_type)_$(AV_type).jld2"

filepath = joinpath(folder, filename)

sol_LDG = solve_ode(u, tspan, params)

dvPdv, maxdvPdv_ind = compute_dvPdv(sol_LDG, equations, md, rd)

common_size = (800, 600)
margins = (5Plots.mm, 5Plots.mm, 5Plots.mm, 5Plots.mm) 

p = plot(sol_LDG.t, dvPdv, 
    xlabel = L"$t$", 
    #ylabel = L"$\max_k{\left(\frac{||\delta v||^2}{||\Pi_N \delta v||^2}\right)}$", 
    label = false , 
    linewidth = 3, tickfontsize = 20,
    legendfontsize = 16,
    xguidefontsize = 24,
    size = common_size,
    left_margin = margins[1],
    right_margin = margins[2],
    top_margin = margins[3],
    bottom_margin = margins[4])

path_to_plots = joinpath(@__DIR__, "Shu_Osher_plots")
mkpath(path_to_plots)
savefig(p, joinpath(path_to_plots, "dvPdv_$(sol_LDG.prob.p.DG_type).pdf"))
# display(p)

AV_type = :BR1
t = Float64[]
max_epsilon = Float64[] 
du_visc_norm = Float64[]
sigma_norm = Float64[]
du_Euler_norm = Float64[]
du_Euler_AV_norm = Float64[]
params = (; rd, md, equations, interface_flux = flux_hllc, 
        initial_condition_type, AV_type, t, max_epsilon, DG_type, 
        du_visc_norm, sigma_norm, du_Euler_norm, du_Euler_AV_norm)
filename = "N$(N)_K$(K1D)_tspan$(tspan[2])_$(DG_type)_$(AV_type).jld2"

filepath = joinpath(folder, filename)

sol_BR1 = solve_ode(u, tspan, params)



du_LDG_compare = plot(sol_LDG.prob.p.t, sol_LDG.prob.p.du_Euler_AV_norm,
    xlabel = L"$t$", ylabel = L"$||du||^2$",
    ylims = (400,1000),
    size = common_size,
    left_margin = margins[1],
    right_margin = margins[2],
    top_margin = margins[3],
    bottom_margin = margins[4],
    title = "$(DG_type) DG with $(sol_LDG.prob.p.AV_type) AV",
    label = "Euler + AV du" ,
    linewidth = 2, guidefontsize = 14,   # for axis labels
    tickfontsize = 12,
    fontsize = 14)

plot!(du_LDG_compare, sol_LDG.prob.p.t, sol_LDG.prob.p.du_Euler_norm,
    label = "Euler + AV du", linewidth = 2) 
    

sigma_norm_evolution_plot = plot(sol_LDG.prob.p.t, sol_LDG.prob.p.sigma_norm,
    xlabel = L"$t$", ylabel = L"$||\sigma||^2$",
    size = common_size,
    left_margin = margins[1],
    right_margin = margins[2],
    top_margin = margins[3],
    bottom_margin = margins[4],
    title = "$(DG_type) DG",
    label = "$(sol_LDG.prob.p.AV_type) $(sol_LDG.destats.naccept) time steps" ,
    linewidth = 2, guidefontsize = 14,   # for axis labels
    tickfontsize = 12,
    fontsize = 14)

plot!(sigma_norm_evolution_plot, sol_BR1.prob.p.t, sol_BR1.prob.p.sigma_norm,
    label = "$(sol_BR1.prob.p.AV_type) $(sol_BR1.destats.naccept) time steps")

savefig(sigma_norm_evolution_plot, joinpath(path_to_plots, "sigma_norm_evolution_$(DG_type).pdf"))

du_visc_norm_evolution_plot = plot(sol_LDG.prob.p.t[1:end-1], diff(sol_LDG.prob.p.du_visc_norm),
    xlabel = L"$t$", ylabel = L"$||du_{visc}||^2$", ylims = (-100.0,100.0),
    size = common_size,
    left_margin = margins[1],
    right_margin = margins[2],
    top_margin = margins[3],
    bottom_margin = margins[4],
    title = "$(DG_type) DG",
    label = "$(sol_LDG.prob.p.AV_type) $(sol_LDG.destats.naccept) time steps" ,
    linewidth = 2, guidefontsize = 14,   # for axis labels
    tickfontsize = 12,
    fontsize = 14)

plot!(du_visc_norm_evolution_plot, sol_BR1.prob.p.t[1:end-1], diff(sol_BR1.prob.p.du_visc_norm),
    label = "$(sol_BR1.prob.p.AV_type) $(sol_BR1.destats.naccept) time steps")

savefig(du_visc_norm_evolution_plot, joinpath(path_to_plots, "du_visc_norm_evolution_$(DG_type).pdf"))



epsilon_evolution_plot = plot(sol_LDG.prob.p.t,abs.(sol_LDG.prob.p.max_epsilon) .+ 1e-14,
    yscale = :log10,
    ylim = (1e-3,1e-1),
    xlabel = L"$t$", 
    #ylabel = L"$\log_{10}\left(\max_k{\epsilon}\right)$",
    size = common_size,
    left_margin = margins[1],
    right_margin = margins[2],
    top_margin = margins[3],
    bottom_margin = margins[4],
    label = "$(sol_LDG.prob.p.AV_type) $(sol_LDG.destats.naccept) time steps" ,
    linewidth = 3, tickfontsize = 20,
    legendfontsize = 16,
    xguidefontsize = 24,)

plot!(epsilon_evolution_plot, sol_BR1.prob.p.t, abs.(sol_BR1.prob.p.max_epsilon) .+ 1e-14,
    linewidth = 3,
    yscale = :log10,
    label = "$(sol_BR1.prob.p.AV_type) $(sol_BR1.destats.naccept) time steps")

savefig(epsilon_evolution_plot, joinpath(path_to_plots, "epsilon_evolution_$(DG_type).pdf"))

up = rd.Vp * getindex.(parent(sol_LDG.u[end]), 1)
xp = rd.Vp * md.x

density_plot = plot(vec(xp), vec(up), label = "$(sol_LDG.prob.p.AV_type) AV", 
            xlabel = L"$x$", 
            #ylabel=L"$\rho$", 
            # ylims = (1e-7,1e0),
            size = common_size,
            left_margin = margins[1],
            right_margin = margins[2],
            top_margin = margins[3],
            bottom_margin = margins[4],
            linewidth = 3, tickfontsize = 20,
            legendfontsize = 16,
            xguidefontsize = 24,
            legend =:topright)

up = rd.Vp * getindex.(parent(sol_BR1.u[end]), 1)

plot!(density_plot, vec(xp), vec(up),
linewidth = 3,
    label = "$(sol_BR1.prob.p.AV_type) AV")

savefig(density_plot, joinpath(path_to_plots, "density_$(sol_LDG.prob.p.DG_type).pdf"))
