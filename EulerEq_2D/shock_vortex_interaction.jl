using StartUpDG
using LinearAlgebra 
using StaticArrays
using OrdinaryDiffEq, RecursiveArrayTools
using JLD2
using Trixi
using Trixi.ForwardDiff
using Plots

# initial_condition_type = :Shock_vortex_interaction

include("initial_conditions.jl")
include("rhs.jl")

psi_1(u, normal, ::CompressibleEulerEquations2D) = u[2] * normal
psi_2(u, normal, ::CompressibleEulerEquations2D) = u[3] * normal

dudv(v, equations) = ForwardDiff.jacobian(v -> entropy2cons(v, equations), v)

gamma = 1.4
equations = CompressibleEulerEquations2D(gamma)

initial_condition_type = :Shock_vortex_interaction

#Choose Modal or Nodal
#DG_type = :nodal
DG_type = :modal

#Mesh Parameters
N = 2 #Order of Polynomial
K1D = 64 #Number of elements in each direction

#BR1 versus LDG switch
#AV_discretization = :BR1
AV_discretization = :LDG

if DG_type == :nodal
        rd = RefElemData(Tri(), SBP(), N)
elseif DG_type == :modal
        rd = RefElemData(Tri(), N)
end

(VX, VY), EToV = uniform_mesh(rd.element_type, 2*K1D, K1D)
VX = VX .+ 1.0
VY = 0.5 * (VY .+ 1.0)
md = MeshData((VX, VY), EToV, rd; is_periodic = (true, false))

u_init = rd.Pq * Shock_vortex_interaction.(md.xq, md.yq, equations, gamma)

Final_time = 0.7
tspan = (0.0, Final_time) 
params = (; rd, md, equations, interface_flux = flux_hllc, AV_discretization, gamma)

ode = ODEProblem(rhs!, VectorOfArray(u_init), tspan, params)
sol = solve(ode, SSPRK43(); abstol = 1e-6, reltol = 1e-4, 
            saveat = LinRange(tspan..., 50), 
            callback = AliveCallback(alive_interval=100))

data_path = joinpath(@__DIR__, "Shock_vortex_interaction_Data")
mkpath(data_path) 
xp, yp = md.x, md.y

# Define the indices we want to plot: Start, Middle-ish, and End
plot_indices = [1, length(sol.t)-35, length(sol.t)]

for i in plot_indices
    u_curr = parent(sol.u[i])
    t_val = sol.t[i]
    rho = getindex.(u_curr, 1)

    # --- 1. Calculate & Plot Schlieren ---
    drhodx = (md.rxJ .* (rd.Dr * rho) + md.sxJ .* (rd.Ds * rho)) ./ md.J
    drhody = (md.ryJ .* (rd.Dr * rho) + md.syJ .* (rd.Ds * rho)) ./ md.J
    g = sqrt.(drhodx.^2 + drhody.^2)
    
    gmax, gmin = maximum(g), minimum(g)
    rho_schl = @. exp(-10 * (g - gmin) / (gmax - gmin))
    plots_path = joinpath(@__DIR__, "Shock_vortex_interaction_plots")
    mkpath(plots_path)
    p_schl = scatter(vec(xp), vec(yp), zcolor=vec(rho_schl),
                        msw=0, ms=0.5, legend=false, ratio=1, cam=(0,90), 
                        c=:red, colorbar=true, guidefontsize=20, tickfontsize=20)
    savefig(p_schl, joinpath(plots_path, "schlieren_t_$(t_val).png"))
    export_to_vtk(rd, md, [rho_schl], joinpath(plots_path, "schlieren_t_$(t_val).vtu"))
    

    # --- 2. Plot Density ---
    p_dens = scatter(vec(xp), vec(yp), zcolor=vec(rho), clims=(0.85, 1.20), 
                        msw=0, ms=0.5, legend=false, ratio=1, cam=(0,90),
                        c=:red, colorbar=true, guidefontsize=20, tickfontsize=20)
    savefig(p_dens, joinpath(plots_path, "density_t_$(t_val).png"))
    export_to_vtk(rd, md, [rho], joinpath(plots_path, "density_t_$(t_val).vtu"))
end

# --- 3. Save Data ---
filename = joinpath(data_path, "N$(N)_K$(2*K1D)x$(K1D)_T$(Final_time)_$(DG_type)_$(AV_discretization).jld2")
@save filename md rd sol
