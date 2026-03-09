using StartUpDG
using LinearAlgebra
using OrdinaryDiffEq, RecursiveArrayTools
using JLD2
using Trixi
using Trixi.ForwardDiff
using StaticArrays
using Plots
using LaTeXStrings

initial_condition_type = :density_wave

include("initial_conditions.jl")
include("rhs.jl")


#Choose AV type
# AV_type = :BR1
AV_type = :LDG

#Choose modal or nodal DG
DG_type = :modal
#DG_type = :nodal

gamma = 1.4
equations = CompressibleEulerEquations1D(gamma)

N = 5
K1D = 16  
if DG_type == :modal
    rd = RefElemData(Line(), N) #modal
elseif DG_type == :nodal
    rd = RefElemData(Line(),SBP(), N) #nodal
end

(VX,), EToV = uniform_mesh(rd.element_type, K1D)



md = MeshData((VX,),EToV, rd; is_periodic=true)
x = md.x
u = rd.Pq * my_initial_condition_density_wave.(md.xq, equations) 

folder = joinpath(@__DIR__, "Density_Wave_data")
    mkpath(folder)
    AV_type = :LDG
    t = Float64[]
    l2_epsilon = Float64[]
    max_epsilon = Float64[]
    du_visc_norm = Float64[]
    sigma_norm = Float64[]
    du_Euler_norm = Float64[]
    du_Euler_AV_norm = Float64[]
    L2_error = Float64[]
    tspan = (0.0, 25.0)
    params = (; rd, md, equations, interface_flux = flux_hllc, 
        initial_condition_type, AV_type, t, max_epsilon, DG_type,
        du_visc_norm, sigma_norm, du_Euler_norm, du_Euler_AV_norm, L2_error, l2_epsilon)
    filename = "N$(N)_K$(K1D)_tspan$(tspan[2])_$(DG_type)_$(AV_type).jld2"
    
    filepath = joinpath(folder, filename)

    ode = ODEProblem(rhs!, VectorOfArray(u), tspan, params)
    sol = solve(ode, SSPRK43(); abstol = 1e-6, reltol = 1e-4, 
                saveat=LinRange(tspan..., 1000), 
                callback=AliveCallback(alive_interval=100))

    t = Float64[]
    l2_epsilon = Float64[]
    max_epsilon = Float64[]
    du_visc_norm = Float64[]
    sigma_norm = Float64[]
    du_Euler_norm = Float64[]
    du_Euler_AV_norm = Float64[]
    L2_error = Float64[]
    params = (; rd, md, equations, interface_flux = flux_hllc, 
        initial_condition_type, AV_type, t, max_epsilon, DG_type,
        du_visc_norm, sigma_norm, du_Euler_norm, du_Euler_AV_norm, L2_error, l2_epsilon)

    ode = ODEProblem(rhs_shock_cap!, VectorOfArray(u), tspan, params)
    sol_sc = solve(ode, SSPRK43(); abstol = 1e-6, reltol = 1e-4, 
                saveat=LinRange(tspan..., 1000), 
                callback=AliveCallback(alive_interval=100))

    t = Float64[]
    du_Euler_norm = Float64[]
    L2_error = Float64[]
    params = (; rd, md, equations, interface_flux = flux_hllc, 
        initial_condition_type, AV_type, t, max_epsilon, DG_type,
        du_visc_norm, sigma_norm, du_Euler_norm, du_Euler_AV_norm, L2_error, l2_epsilon)

    ode = ODEProblem(rhs_DG_FEM!, VectorOfArray(u), tspan, params)
    sol_DG_FEM = solve(ode, SSPRK43(); abstol = 1e-6, reltol = 1e-4, 
                saveat=LinRange(tspan..., 1000), 
                callback=AliveCallback(alive_interval=100))

    common_size = (800, 600)
    margins = (5Plots.mm, 5Plots.mm, 5Plots.mm, 5Plots.mm)   # (left, right, top, bottom)
    
    #Plot the maximum epsilon values
    p = plot(sol.prob.p.t, sol.prob.p.max_epsilon .+ 1e-14, 
        yscale = :log10,
        ylims = (1e-11,1e0),
        size = common_size,
        yticks = [1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e0],
        left_margin = margins[1],
        right_margin = margins[2],
        top_margin = margins[3],
        bottom_margin = margins[4],
        # title = "Density wave N = $(rd.N) K = $(md.K)",
        xlabel = L"$t$",
        linewidth = 3, tickfontsize = 20,
        legendfontsize = 16,
        xguidefontsize = 24,
        legend=:topright,
        label = L"$\max_k{\epsilon_k}$ ECAV")

    plot!(p, sol_sc.prob.p.t, sol_sc.prob.p.max_epsilon .+ 1e-14, 
        label = L"$\max_k{\epsilon_k}$ SC", linewidth = 3)

    display(p)

    path_to_plots = joinpath(@__DIR__, "Density_Wave_Plots")
    mkpath(path_to_plots)
    
    savefig(p, joinpath(path_to_plots, "max_epsilon_$(sol.prob.p.DG_type)_K_$(K1D).png"))

    #Plot the l2 norm of the epsilon values

    p = plot(sol.prob.p.t, sol.prob.p.l2_epsilon .+ 1e-14, 
        yscale = :log10,
        ylims = (1e-15,1e0),
        size = common_size,
        left_margin = margins[1],
        right_margin = margins[2],
        top_margin = margins[3],
        bottom_margin = margins[4],
        # title = "Density wave N = $(rd.N) K = $(md.K)",
        xlabel = L"$t$",
        linewidth = 3, tickfontsize = 20,
        legendfontsize = 16,
        xguidefontsize = 24,
        legend=:topright,
        label = label = L"$\||\epsilon\||_{\ell^2}$ ECAV")

    plot!(p, sol_sc.prob.p.t, sol_sc.prob.p.l2_epsilon .+ 1e-14, 
        label = label = L"$\||\epsilon\||_{\ell^2}$ SC", linewidth = 3)

    display(p)
    savefig(p, joinpath(path_to_plots, "l2_epsilon_$(sol.prob.p.DG_type)_K_$(K1D).png"))

    #Plot the L2 error evolution
    q = plot(sol.prob.p.t, sol.prob.p.L2_error,
        # title = "Density wave N = $(rd.N) K = $(md.K)",
        xlabel = L"$t$",
        size = common_size,
        left_margin = margins[1],
        right_margin = margins[2],
        top_margin = margins[3],
        bottom_margin = margins[4],
        label = L"$L^2$ error ECAV", linewidth = 3,
        tickfontsize = 20,
        legendfontsize = 16,
        xguidefontsize = 24,
        legend=:right
    )

    plot!(q,sol_sc.prob.p.t, sol_sc.prob.p.L2_error,
        label = L"$L^2$ error SC", linewidth = 3)

    display(q)

    u_ecav = sol.u[end]
    u_DG_FEM = sol_DG_FEM.u[end]
    error_ecav_dg_fem = rd.Vq * (parent(u_ecav) - parent(u_DG_FEM))
    L2_error_ecav_dg_fem = sqrt(sum(md.wJq .* norm.(error_ecav_dg_fem).^2))

    print("The L2 norm of the difference between "*
     "the ECAV DG and DG FEM solution at the final "*
     "time is: $(L2_error_ecav_dg_fem)\n "*
     "The error for the ECAV solution at the final time is: $(sol.prob.p.L2_error[end]) \n "*
     "The error for the DG FEM solution at the final time is: $(sol_DG_FEM.prob.p.L2_error[end]).\n ")

    # Compare the final time errors in a latex table format
    sig = 4
    latex_row = "$(K1D) & " * "\$$(round(sol.prob.p.L2_error[end], sigdigits = sig))\$ & " *
    "\$$(round(sol_DG_FEM.prob.p.L2_error[end], sigdigits = sig))\$ & " *
    "\$$(round(L2_error_ecav_dg_fem, sigdigits = sig))\$ \\\\"
    println(latex_row)

    savefig(q, joinpath(path_to_plots, "error_$(sol.prob.p.DG_type)_K_$(K1D).png"))

    #Plot densities at the final time
    rho_sc = rd.Vp * getindex.(parent(sol_sc.u[end]), 1)
    rho_ecav = rd.Vp * getindex.(parent(sol.u[end]), 1)
    x_refined = range(minimum(vec(md.x)), maximum(vec(md.x)), 10*length(vec(md.x)))
    u_exact = my_initial_condition_density_wave.(x_refined .- 0.1*tspan[2], equations)
    rho_exact = getindex.(u_exact, 1)

    x_plot = vec(rd.Vp * md.x)

    p = plot(x_plot, vec(rho_ecav), 
                linewidth = 3,
                label = "",
                color = :blue,
                size = common_size,
                left_margin = margins[1],
                right_margin = margins[2],
                top_margin = margins[3],
                bottom_margin = margins[4],
                tickfontsize = 20,
                legendfontsize = 16,
                xguidefontsize = 24,
                xlabel = L"$x$",
                legend = :topright
                )

    inc = 1

    sx = x_plot[1:inc:end]
    sy = vec(rho_ecav)[1:inc:end]

    scatter!(p,sx,sy, marker = :circle, ms = 4, color =:blue, label = "", primary = false)

    plot!(p, x_plot, vec(rho_sc), linewidth = 3,   
        label = "")

    inc = 5
    sx = x_plot[1:inc:end]
    sy = vec(rho_sc)[1:inc:end]

    scatter!(p,sx,sy, marker = :rect, ms = 4, color =:red, label = "", primary = false)
    
    plot!(
    p,
    Float64[], Float64[],
    linewidth = 3,
    marker = :circle,
    ms = 4,
    color = :blue,
    label = "ECAV"
    )

    plot!(
    p,
    Float64[], Float64[],
    linewidth = 3,
    marker = :rect,
    ms = 4,
    color = :red,
    label = "SC"
    )

    plot!(p, vec(collect(x_refined)), vec(rho_exact), linewidth = 3,
        label = "Exact",
        color = :green)
    display(p)

    savefig(p, joinpath(path_to_plots, "density_$(sol.prob.p.DG_type)_K_$(K1D).png"))


elseif initial_condition_type == :sod_shock_tube
    
    filepath = joinpath(folder, filename)

    if isfile(filepath)
        data = load(filepath)
        md = data["md"]
        rd = data["rd"]
        sol = data["sol"]
    else
        sol = solve_ode(u, tspan, params)
        @save filepath md rd sol
    end

    dvPdv = compute_dvPdv(sol, equations, md, rd)

    p = plot(sol.t, dvPdv, 
        title = "Modified Sod shock tube",
        xlabel = L"$t$", ylabel = L"$\max_k{\left(\frac{||\delta v||^2}{||\Pi_N \delta v||^2}\right)}$" 
        , left_margin = 10Plots.mm,
        label = "N = $(rd.N) K = $(md.K) " , linewidth = 2)

    display(p)