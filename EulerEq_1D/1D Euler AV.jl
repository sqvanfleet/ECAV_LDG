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
include("computing_quanities.jl")

psi(u, normal, ::CompressibleEulerEquations1D) = u[2] * normal
dudv(v, equations) = ForwardDiff.jacobian(v -> entropy2cons(v, equations), v)

function ramp(s,s0,xi,epsilon_0)
    if s < s0 - xi
        return 0.0
    elseif abs(s-s0) < xi
        return epsilon_0/2 * (1 + sin(pi/(2*xi) * (s - s0)))
    else
        return epsilon_0
    end
end

function rhs_shock_cap!(du_voa, u_voa, params, t)
    du = parent(du_voa)
    u = parent(u_voa)

    (; rd, md, equations) = params

    push!(params.t, t)

    #compute entropy variables 
    v = rd.Pq * cons2entropy.(rd.Vq * u, equations)

    invMQTr = -rd.M \ (rd.Dr' * rd.M)
    # uM = rd.Vf * u # nodal
    uM = entropy2cons.(rd.Vf * v, equations) # modal 
    uP = uM[md.mapP]
    if params.initial_condition_type == :sod_shock_tube
        uP[1] = initial_condition_sod_shock_tube(md.x[1], equations)
        uP[end] = initial_condition_sod_shock_tube(md.x[end], equations)
    elseif params.initial_condition_type == :shu_osher
        uP[1] = shu_osher(md.x[1], equations)
        uP[end] = shu_osher(md.x[end], equations)
    end
    # vM = rd.Vf * v #nodal
    # vP = vM[md.mapP] #nodal
    vM = cons2entropy.(uM, equations) #modal
    vP = cons2entropy.(uP, equations) #modal
    interface_flux = @. params.interface_flux(uM, uP, SVector(md.nxJ), equations) * md.Jf
    duvol = md.rxJ .* (invMQTr * rd.Pq * flux.(rd.Vq * u, 1, equations)) 
    du .= duvol + rd.LIFT * interface_flux

    du_Euler_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * (du./md.J)).^2))
    push!(params.du_Euler_norm, du_Euler_norm)

    # use ρ as indicator
    rho = getindex.(u, 1)
    indicator = pressure.(u, equations) .* rho
    u_modal = rd.VDM \ indicator
    EN = vec(sum(u_modal.^2, dims=1))
    ENm1 = vec(sum(u_modal[1:end-1, :].^2, dims=1))
    S = max.(u_modal[end, :].^2 ./ EN, u_modal[end-1, :].^2 ./ ENm1)
    s = log10.(S);
    element_sizes = transpose(maximum(md.x, dims=1) .- minimum(md.x, dims=1))
    epsilon_0 = 0.5*element_sizes/rd.N
    s_0 = -(15/2)*log10(rd.N)
    xi = (7/2)*log10(rd.N)

    epsilon = @. ramp(s, s_0, xi, epsilon_0)

    if params.AV_type == :BR1
        beta = 0.0
    elseif params.AV_type == :LDG
        beta = @. sign(md.nxJ)
    end


    theta = (md.rxJ .* (invMQTr * v) + rd.LIFT * 
            (@. 0.5 * (vP + vM) * md.nxJ - 0.5 * beta * (vP - vM) * md.nxJ)) ./ md.J
    
    K = dudv.(cons2entropy.(rd.Vq * u, equations), equations)

    sigma = rd.Pq * (K .* (rd.Vq * theta))

    sigma = sigma * Diagonal(vec(epsilon))
    sigma_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * sigma).^2))
    
    sigmaM = rd.Vf * sigma
    sigmaP = sigmaM[md.mapP]

    du_visc = md.rxJ .* (invMQTr * sigma) + 
        rd.LIFT * (@. 0.5 * (sigmaP + sigmaM) * md.nxJ + 0.5 * beta * (sigmaP - sigmaM) * md.nxJ) ./ md.J
    
    du_visc_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * du_visc).^2))
    
    du .-= md.rxJ .* (invMQTr * sigma) + 
        rd.LIFT * (@. 0.5 * (sigmaP + sigmaM) * md.nxJ + 0.5 * beta * (sigmaP - sigmaM) * md.nxJ)
    @. du /= -md.J

    du_Euler_AV_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * du).^2))

    if params.initial_condition_type == :density_wave
        u_exact = initial_condition_density_wave.(md.xq .- 0.1*t, equations)
        error = rd.Vq * parent(u) - u_exact
        num = sqrt(sum(md.wJq .* norm.(error).^2))
        den = sqrt(sum(md.wJq .* norm.(u_exact).^2))
        L2_error = num
        push!(params.L2_error, L2_error)
        l2_epsilon = sqrt(sum(epsilon.^2))
        push!(params.l2_epsilon, l2_epsilon)
    end

    push!(params.du_Euler_AV_norm, du_Euler_AV_norm)
    push!(params.max_epsilon, maximum(epsilon))
    push!(params.du_visc_norm, du_visc_norm)
    push!(params.sigma_norm, sigma_norm)

    return du

end

function rhs!(du_voa, u_voa, params, t)
    
    du = parent(du_voa)
    u = parent(u_voa)

    (; rd, md, equations) = params

    push!(params.t, t)

    #compute entropy variables 
    v = rd.Pq * cons2entropy.(rd.Vq * u, equations)

    invMQTr = -rd.M \ (rd.Dr' * rd.M)
    # uM = rd.Vf * u # nodal
    uM = entropy2cons.(rd.Vf * v, equations) # modal 
    uP = uM[md.mapP]
    if params.initial_condition_type == :sod_shock_tube
        uP[1] = initial_condition_sod_shock_tube(md.x[1], equations)
        uP[end] = initial_condition_sod_shock_tube(md.x[end], equations)
    elseif params.initial_condition_type == :shu_osher
        uP[1] = shu_osher(md.x[1], equations)
        uP[end] = shu_osher(md.x[end], equations)
    end
    # vM = rd.Vf * v #nodal
    # vP = vM[md.mapP] #nodal
    vM = cons2entropy.(uM, equations) #modal
    vP = cons2entropy.(uP, equations) #modal
    interface_flux = @. params.interface_flux(uM, uP, SVector(md.nx), equations) * md.Jf
    duvol = md.rxJ .* (invMQTr * rd.Pq * flux.(rd.Vq * u, 1, equations)) 
    du .= duvol + rd.LIFT * interface_flux

    du_Euler_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * (du./md.J)).^2))
    push!(params.du_Euler_norm, du_Euler_norm)

    if params.AV_type == :BR1
        beta = 0.0
    elseif params.AV_type == :LDG
        beta = @. sign(md.nxJ)
    end

    theta = (md.rxJ .* (invMQTr * v) + rd.LIFT * 
            (@. 0.5 * (vP + vM) * md.nxJ - 0.5 * beta * (vP - vM) * md.nxJ)) ./ md.J
    term1 = dot.(v, rd.M * duvol)
    term2 = @. rd.wf * psi(uM, md.nxJ, equations)
    delta = sum(term1, dims=1) + sum(term2, dims=1)
    num = @. -min(0, delta)
    v_avg = repeat(0.5 * sum(Diagonal(rd.wq) * (rd.Vq * v), dims=1), rd.Np, 1)
    K = dudv.(cons2entropy.(rd.Vq * u, equations), equations)
    #K = dudv.(rd.Vq * v_avg, equations)
    sigma = rd.Pq * (K .* (rd.Vq * theta))
    den = sum(md.wJq .* dot.(rd.Vq * sigma, rd.Vq * theta), dims=1)
    epsilon = @. num * den / (1e-14 + den^2)
    sigma = sigma * Diagonal(vec(epsilon))
    sigma_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * sigma).^2))
    
    sigmaM = rd.Vf * sigma
    sigmaP = sigmaM[md.mapP]

    du_visc = md.rxJ .* (invMQTr * sigma) + 
        rd.LIFT * (@. 0.5 * (sigmaP + sigmaM) * md.nxJ + 0.5 * beta * (sigmaP - sigmaM) * md.nxJ) ./ md.J
    
    du_visc_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * du_visc).^2))
    
    du .-= md.rxJ .* (invMQTr * sigma) + 
        rd.LIFT * (@. 0.5 * (sigmaP + sigmaM) * md.nxJ + 0.5 * beta * (sigmaP - sigmaM) * md.nxJ)
    @. du /= -md.J

    du_Euler_AV_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * du).^2))

    if params.initial_condition_type == :density_wave
        u_exact = initial_condition_density_wave.(md.xq .- 0.1*t, equations)
        error = rd.Vq * parent(u) - u_exact
        num = sqrt(sum(md.wJq .* norm.(error).^2))
        den = sqrt(sum(md.wJq .* norm.(u_exact).^2))
        L2_error = num
        push!(params.L2_error, L2_error)
        l2_epsilon = sqrt(sum(epsilon.^2))
        push!(params.l2_epsilon, l2_epsilon)
    end

    if params.initial_condition_type == :density_wave_convergence
        u_exact = density_wave_convergence.(md.xq .- 0.1*t, equations)
        error = rd.Vq * parent(u) - u_exact
        num = sqrt(sum(md.wJq .* norm.(error).^2))
        den = sqrt(sum(md.wJq .* norm.(u_exact).^2))
        L2_error = num/den
        push!(params.L2_error, L2_error)
    end

    push!(params.du_Euler_AV_norm, du_Euler_AV_norm)
    push!(params.max_epsilon, maximum(epsilon))
    push!(params.du_visc_norm, du_visc_norm)
    push!(params.sigma_norm, sigma_norm)

    return du
end

function rhs_DG_FEM!(du_voa, u_voa, params, t)
    
    du = parent(du_voa)
    u = parent(u_voa)

    (; rd, md, equations) = params

    push!(params.t, t)

    #compute entropy variables 
    v = rd.Pq * cons2entropy.(rd.Vq * u, equations)

    invMQTr = -rd.M \ (rd.Dr' * rd.M)
    # uM = rd.Vf * u # nodal
    uM = entropy2cons.(rd.Vf * v, equations) # modal 
    uP = uM[md.mapP]
    if params.initial_condition_type == :sod_shock_tube
        uP[1] = initial_condition_sod_shock_tube(md.x[1], equations)
        uP[end] = initial_condition_sod_shock_tube(md.x[end], equations)
    elseif params.initial_condition_type == :shu_osher
        uP[1] = shu_osher(md.x[1], equations)
        uP[end] = shu_osher(md.x[end], equations)
    end

    vM = cons2entropy.(uM, equations) 
    vP = cons2entropy.(uP, equations) 
    interface_flux = @. params.interface_flux(uM, uP, SVector(md.nx), equations) * md.Jf
    duvol = md.rxJ .* (invMQTr * rd.Pq * flux.(rd.Vq * u, 1, equations)) 
    du .= duvol + rd.LIFT * interface_flux

    du_Euler_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * (du./md.J)).^2))
    push!(params.du_Euler_norm, du_Euler_norm)

    @. du /= -md.J

    if params.initial_condition_type == :density_wave
        u_exact = initial_condition_density_wave.(md.xq .- 0.1*t, equations)
        error = rd.Vq * parent(u) - u_exact
        num = sqrt(sum(md.wJq .* norm.(error).^2))
        den = sqrt(sum(md.wJq .* norm.(u_exact).^2))
        L2_error = num
        push!(params.L2_error, L2_error)
    end

    if params.initial_condition_type == :density_wave_convergence
        u_exact = density_wave_convergence.(md.xq .- 0.1*t, equations)
        error = rd.Vq * parent(u) - u_exact
        num = sqrt(sum(md.wJq .* norm.(error).^2))
        den = sqrt(sum(md.wJq .* norm.(u_exact).^2))
        L2_error = num/den
        push!(params.L2_error, L2_error)
    end

    return du
end

#Choose AV type
# AV_type = :BR1
AV_type = :LDG

#Choose modal or nodal DG
DG_type = :modal
#DG_type = :nodal

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

#Choose the problem
# initial_condition_type = :density_wave
# initial_condition_type = :density_wave_convergence
# initial_condition_type = :sod_shock_tube
# initial_condition_type = :stationary_contact_wave
initial_condition_type = :shu_osher
# initial_condition_type = :Euler_problem_1

if initial_condition_type == :density_wave
    md = MeshData((VX,),EToV, rd; is_periodic=true)
    x = md.x
    u = rd.Pq * initial_condition_density_wave.(md.xq, equations) 
    
elseif initial_condition_type == :sod_shock_tube
    @. VX = 0.5 * (VX .+ 1);
    md = MeshData((VX,), EToV, rd)
    x = md.x
    u = rd.Pq * initial_condition_sod_shock_tube.(md.xq, equations)

elseif initial_condition_type == :stationary_contact_wave
    md = MeshData((VX,), EToV, rd; is_periodic=true)
    x = md.x
    u = rd.Pq * stationary_contact_wave.(md.xq, equations)

elseif initial_condition_type == :shu_osher
    @. VX = 5.0*VX
    md = MeshData((VX,), EToV, rd)
    u = rd.Pq * shu_osher.(md.xq, equations)

end

# params = (; rd, md, equations, interface_flux = flux_hllc, initial_condition_type, AV_type)

function solve_ode(u, tspan, params, abstol = 1e-6, reltol = 1e-4, number_of_saves = 1000)
    # Solve the ODE problem
    ode = ODEProblem(rhs!, VectorOfArray(u), tspan, params)
    sol = solve(ode, SSPRK43(); abstol, reltol, 
                saveat=LinRange(tspan..., number_of_saves), 
                callback=AliveCallback(alive_interval=100))
    return sol
end

if initial_condition_type == :density_wave_convergence
    K_list = [8,16,32,64]
    N_list = [1,2,3,4]

    for N in N_list
        for K1D in K_list
            if DG_type == :modal
                local rd = RefElemData(Line(), N) #modal
            elseif DG_type == :nodal
                local rd = RefElemData(Line(),SBP(), N) #nodal
            end
            local (VX,), EToV = uniform_mesh(rd.element_type, K1D)
            local md = MeshData((VX,),EToV, rd; is_periodic=true)
            local u = rd.Pq * density_wave_convergence.(md.xq, equations) 
            local tspan = (0.0, 10.0)
            local folder = "/Users/samvanfleet/Documents/Rice Artificial Viscosity/Euler_1D copy/Density_Wave_convergence_data/"
            local AV_type = :LDG
            local t = Float64[]
            local max_epsilon = Float64[]
            local du_visc_norm = Float64[]
            local sigma_norm = Float64[]
            local du_Euler_norm = Float64[]
            local du_Euler_AV_norm = Float64[]
            local L2_error = Float64[]
            local params = (; rd, md, equations, interface_flux = flux_hllc, 
                initial_condition_type, AV_type, t, max_epsilon, DG_type,
                du_visc_norm, sigma_norm, du_Euler_norm, du_Euler_AV_norm, L2_error)
            local filename = "N$(N)_K$(K1D)_tspan$(tspan[2])_$(DG_type)_$(AV_type).jld2"

            local filepath = joinpath(folder, filename)

            local sol = solve_ode(u, tspan, params)

            @save filepath md rd sol
        end
    end

elseif  initial_condition_type == :density_wave
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
    
    savefig(p, joinpath("Density_Wave_Plots", "max_epsilon_$(sol.prob.p.DG_type)_K_$(K1D).png"))

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
    savefig(p, joinpath("Density_Wave_Plots", "l2_epsilon_$(sol.prob.p.DG_type)_K_$(K1D).png"))

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

    savefig(q, joinpath("Density_Wave_Plots", "error_$(sol.prob.p.DG_type)_K_$(K1D).png"))

    #Plot densities at the final time
    rho_sc = rd.Vp * getindex.(parent(sol_sc.u[end]), 1)
    rho_ecav = rd.Vp * getindex.(parent(sol.u[end]), 1)
    x_refined = range(minimum(vec(md.x)), maximum(vec(md.x)), 10*length(vec(md.x)))
    u_exact = initial_condition_density_wave.(x_refined .- 0.1*tspan[2], equations)
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

    savefig(p, joinpath("Density_Wave_Plots", "density_$(sol.prob.p.DG_type)_K_$(K1D).png"))


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


elseif initial_condition_type == :stationary_contact_wave
    AV_type = :LDG
    tspan = (0.0, 4.0)
    t = Float64[]
    max_epsilon = Float64[]
    du_visc_norm = Float64[]
    sigma_norm = Float64[]
    du_Euler_norm = Float64[]
    du_Euler_AV_norm = Float64[]
    L2_error = Float64[]
    params = (; rd, md, equations, interface_flux = flux_hllc, 
        initial_condition_type, AV_type, t, max_epsilon, DG_type,
        du_visc_norm, sigma_norm, du_Euler_norm, du_Euler_AV_norm, L2_error)

    folder = joinpath(@__DIR__, "Density_Wave_data")

    mkpath(folder)

    filename = "N$(N)_K$(K1D)_tspan$(tspan[2]).jld2"
    
    filepath = joinpath(folder, filename)

    ode = ODEProblem(rhs!, VectorOfArray(u), tspan, params)
    sol = solve(ode, SSPRK43(); abstol = 1e-6, reltol = 1e-4, 
                saveat=LinRange(tspan..., 1000), 
                callback=AliveCallback(alive_interval=100),
                adaptive = false, 
                dt = 5e-4)
    
    L2_error = compute_L2_error_evolution_stationary_contact_wave(sol, equations, md, rd)

    @save filepath md rd sol L2_error


elseif initial_condition_type == :shu_osher

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
    

elseif initial_condition_type == :Euler_problem_1
    md = MeshData((VX,), EToV, rd)
    pmin_vals = range(1e-3, stop=0.5, length = 100)

    k = 1.0

    dvPdv = zeros(length(pmin_vals))

    for (i,pmin) in enumerate(pmin_vals)
        local u = rd.Pq * Euler_problem_1.(md.xq, equations, k, pmin, pmin)
        local v = cons2entropy.(rd.Vq * parent(u),equations)
        local v_avg = repeat(0.5 * sum(Diagonal(rd.wq) * (rd.Vq * rd.Pq * v), dims=1), rd.Np, 1)
        local num = sum(md.wJq .* norm.(v - rd.Vq * v_avg).^2, dims = 1)
        local den = sum(md.wJq .* norm.(rd.Vq * rd.Pq * (v - rd.Vq * v_avg)).^2, dims = 1)
        global dvPdv[i] = maximum(num./den)
    end

    p = plot(pmin_vals, dvPdv,xlabel = "pmin and rhomin", ylabel = L"$\max_k{\left(\frac{||\delta v||^2}{||\Pi_N \delta v||^2}\right)}$",
    title = "k = 1", legend = false, linewidth = 2)

    display(p)

end






