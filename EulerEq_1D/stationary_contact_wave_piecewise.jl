using StartUpDG
using LinearAlgebra
using OrdinaryDiffEq, RecursiveArrayTools
using JLD2
using Trixi
using Trixi.ForwardDiff
using StaticArrays
using Plots
using LaTeXStrings

include("initial_conditions.jl")
include("rhs.jl")
include("computing_quanities.jl")

initial_condition_type = :stationary_contact_wave

function error_study()
    
    AV_type = :LDG 
    DG_type = :modal 
    gamma = 1.4
    equations = CompressibleEulerEquations1D(gamma)
    K1D = 80
    tspan = (0.0, 4.0)

    folder = joinpath(@__DIR__, "Contact_wave_data")
    mkpath(folder)

    # --- Loop over N values ---
    for N in [2, 4, 6, 8]
        println("\n" * "-"^20)
        println("Starting N = $N")
        
        if DG_type == :modal
            rd = RefElemData(Line(), N)
        elseif DG_type == :nodal
            rd = RefElemData(Line(), SBP(), N)
        end

        (VX,), EToV = uniform_mesh(rd.element_type, K1D)
        md = MeshData((VX,), EToV, rd; is_periodic=true)
        
        u = rd.Pq * stationary_contact_wave_piecewise.(md.xq, equations)

        # Re-initialize logging containers for each N
        params = (; rd, md, equations, interface_flux = flux_hllc, 
            initial_condition_type, AV_type, DG_type,
            t = Float64[], 
            max_epsilon = Float64[], 
            du_visc_norm = Float64[], 
            sigma_norm = Float64[], 
            du_Euler_norm = Float64[], 
            du_Euler_AV_norm = Float64[])

        ode = ODEProblem(rhs!, VectorOfArray(u), tspan, params)
        
        sol = solve(ode, SSPRK43(); 
                    abstol = 1e-6, 
                    reltol = 1e-4, 
                    saveat = LinRange(tspan..., 1000), 
                    callback = AliveCallback(alive_interval=200),
                    adaptive = false, 
                    dt = 5e-4)

        L2_error = compute_L2_error_evolution_stationary_contact_wave_piecewise(sol, equations, md, rd)

        filename = "N$(N)_K$(K1D)_tspan$(tspan[2]).jld2"
        filepath = joinpath(folder, filename)
        
        @save filepath md rd sol L2_error
        println("Finished N = $N. Data saved.")
    end
end

error_study()


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

savefig(p1, joinpath(base_path, "Contact_wave_plots", "Contact_wave_error_evolution_piecewise.png"))