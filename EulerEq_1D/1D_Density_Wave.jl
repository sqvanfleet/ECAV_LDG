using StartUpDG
using LinearAlgebra
using OrdinaryDiffEq, RecursiveArrayTools
using JLD2
using Trixi
using Trixi.ForwardDiff
using StaticArrays
using Plots
using LaTeXStrings
using Printf

# --- Setup and Includes ---
initial_condition_type = :density_wave
include("initial_conditions.jl")
include("rhs.jl")

gamma = 1.4
equations = CompressibleEulerEquations1D(gamma)
N = 5
tspan = (0.0, 25.0)
DG_type = :modal 
AV_type = :LDG

# Create output directories
path_to_plots = joinpath(@__DIR__, "Density_Wave_Plots")
mkpath(path_to_plots)

function create_params(rd, md, equations, AV_type, DG_type)
    return (; rd, md, equations, interface_flux = flux_hllc, 
        initial_condition_type, AV_type, DG_type,
        t = Float64[], l2_epsilon = Float64[], max_epsilon = Float64[], 
        du_visc_norm = Float64[], sigma_norm = Float64[], 
        du_Euler_norm = Float64[], du_Euler_AV_norm = Float64[], L2_error = Float64[])
end

# --- Storage for the final table ---
latex_table_rows = String[]

# --- Main Loop over K1D ---
for K1D in [8, 16, 24]
    println("\n" * "="^60)
    println("STARTING SIMULATION: K1D = $K1D")
    println("="^60)

    # DISAMBIGUATION: Use 'local' to prevent global/local scope warnings
    local VX, EToV, rd, md, u_init, sol, sol_sc, sol_DG_FEM
    local common_size, rho_sc, rho_ecav, x_refined, rho_exact, x_plot

    if DG_type == :modal
        rd = RefElemData(Line(), N)
    elseif DG_type == :nodal
        rd = RefElemData(Line(), SBP(), N)
    end

    (VX,), EToV = uniform_mesh(rd.element_type, K1D)
    md = MeshData((VX,), EToV, rd; is_periodic=true)
    u_init = rd.Pq * my_initial_condition_density_wave.(md.xq, equations) 

    # Solver Runs
    sol = solve(ODEProblem(rhs!, VectorOfArray(u_init), tspan, create_params(rd, md, equations, AV_type, DG_type)), 
                SSPRK43(); abstol = 1e-6, reltol = 1e-4, saveat=LinRange(tspan..., 1000), 
                callback=AliveCallback(alive_interval=500))

    sol_sc = solve(ODEProblem(rhs_shock_cap!, VectorOfArray(u_init), tspan, create_params(rd, md, equations, AV_type, DG_type)), 
                   SSPRK43(); abstol = 1e-6, reltol = 1e-4, saveat=LinRange(tspan..., 1000), 
                   callback=AliveCallback(alive_interval=500))

    sol_DG_FEM = solve(ODEProblem(rhs_DG_FEM!, VectorOfArray(u_init), tspan, create_params(rd, md, equations, AV_type, DG_type)), 
                       SSPRK43(); abstol = 1e-6, reltol = 1e-4, saveat=LinRange(tspan..., 1000), 
                       callback=AliveCallback(alive_interval=500))

    # --- Plotting ---
    common_size = (800, 600)
    m_left, m_right, m_top, m_bottom = 12Plots.mm, 5Plots.mm, 5Plots.mm, 12Plots.mm

    # Max Epsilon Plot
    p_eps = plot(sol.prob.p.t, sol.prob.p.max_epsilon .+ 1e-14, 
        yscale = :log10, ylims = (1e-11,1e0), size = common_size,
        left_margin=m_left, bottom_margin=m_bottom,
        yticks = [1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e0],
        xlabel = L"$t$", linewidth = 3, tickfontsize = 16, legendfontsize = 12, xguidefontsize = 20,
        legend=:topright, label = L"$\max_k{\epsilon_k}$ ECAV")
    plot!(p_eps, sol_sc.prob.p.t, sol_sc.prob.p.max_epsilon .+ 1e-14, label = L"$\max_k{\epsilon_k}$ SC", linewidth = 3)
    savefig(p_eps, joinpath(path_to_plots, "max_epsilon_$(DG_type)_K_$(K1D).png"))

    # Plot Final Density Profiles
    rho_sc = rd.Vp * getindex.(parent(sol_sc.u[end]), 1)
    rho_ecav = rd.Vp * getindex.(parent(sol.u[end]), 1)
    x_refined = range(minimum(vec(md.x)), maximum(vec(md.x)), length=500)
    rho_exact = getindex.(my_initial_condition_density_wave.(x_refined .- 0.1*tspan[2], equations), 1)
    x_plot = vec(rd.Vp * md.x)

    p_dens = plot(x_plot, vec(rho_ecav), linewidth = 3, color = :blue, size = common_size,
                ylims = (0.975, 1.525),
                 left_margin=m_left, bottom_margin=m_bottom, tickfontsize = 16, xlabel = L"$x$", legend = :topright, label="")
    scatter!(p_dens, x_plot, vec(rho_ecav), marker=:circle, ms=4, color=:blue, label="", primary=false)
    plot!(p_dens, x_plot, vec(rho_sc), linewidth=3, color=:red, label="")
    scatter!(p_dens, x_plot[1:5:end], vec(rho_sc)[1:5:end], marker=:rect, ms=4, color=:red, label="", primary=false)
    
    # Legend Overlays
    plot!(p_dens, [0], [0], linewidth=3, marker=:circle, color=:blue, label="ECAV")
    plot!(p_dens, [0], [0], linewidth=3, marker=:rect, color=:red, label="SC")
    plot!(p_dens, vec(collect(x_refined)), vec(rho_exact), linewidth=3, color=:green, label="Exact")
    savefig(p_dens, joinpath(path_to_plots, "density_$(DG_type)_K_$(K1D).png"))

    # --- Data Collection for Table ---
    local error_diff, L2_diff
    error_diff = rd.Vq * (parent(sol.u[end]) - parent(sol_DG_FEM.u[end]))
    L2_diff = sqrt(sum(md.wJq .* norm.(error_diff).^2))

    push!(latex_table_rows, @sprintf("%d & \$%.4e\$ & \$%.4e\$ & \$%.4e\$ \\\\", 
          K1D, sol.prob.p.L2_error[end], sol_DG_FEM.prob.p.L2_error[end], L2_diff))
end

# --- Print Final Table ---
println("\n" * "*"^50)
println("FINAL LATEX TABLE")
println("*"^50)
println("\\begin{table}[h!]\n  \\centering\n  \\begin{tabular}{l|ccc}\n    \\hline\n    \$K_{1D}\$ & ECAV \$L^2\$ Error & DG FEM \$L^2\$ Error & \$L^2\$ Diff \\\\\n    \\hline")
for row in latex_table_rows; println("    " * row); end
println("    \\hline\n  \\end{tabular}\n\\end{table}")