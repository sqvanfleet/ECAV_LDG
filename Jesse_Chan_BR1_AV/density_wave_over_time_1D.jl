using StartUpDG
using Plots, LinearAlgebra
using Trixi
using StaticArrays
using OrdinaryDiffEq


equations = CompressibleEulerEquations1D(1.4)

psi(u, ::CompressibleEulerEquations1D) = u[2]
psi(u, normal, ::CompressibleEulerEquations1D) = u[2] * normal

regularized_ratio(a, b; tol=1e-14) = a * b / (b^2 + tol)

using Trixi.ForwardDiff
dudv(v, equations) = ForwardDiff.jacobian(v -> entropy2cons(v, equations), v)

function calc_dg_derivative!(dudx, u, params)
    (; rd, md) = params
    uM = rd.Vf * u
    uP = uM[md.mapP]
    interface_flux = @. 0.5 * (uP - uM) * md.nx
    dudx .= md.rxJ .* (rd.Dr * u) + rd.LIFT * interface_flux
    dudx ./= md.J
end

function rhs!(du, u, params, t)

    (; rd, md, equations, interface_flux) = params

    # compute entropy variables
    (; viscosity) = params
    v = rd.Pq * cons2entropy.(rd.Vq * u, equations)

    # calculate derivatives for viscous terms
    (; sigma) = params 
    calc_dg_derivative!(sigma, v, params)

    u_q = rd.Vq * u
    # u_q = entropy2cons.(rd.Vq * v, equations)
    u_avg = repeat(0.5 * sum(rd.M, dims=1) * u, length(rd.rq), 1)

    if params.use_entropy_projection == true
        uM = entropy2cons.(rd.Vf * v, equations)
    else
        uM = rd.Vf * u
    end
    uP = uM[md.mapP]

    # surface terms
    E = rd.Vf   
    du .= E' * (@. interface_flux(uM, uP, SVector(md.nx), equations))

    # volume terms
    Q_skew = 0.5 * (rd.M * rd.Dr - rd.Dr' * rd.M)
    QTr = rd.Dr' * rd.M * rd.Pq
    du_vol = similar(du[:, 1])
    integrand = zeros(size(u_q, 1))
    for e in axes(du, 2)
        if rd.approximation_type isa SBP && params.use_EC_volume_flux == true
            fill!(du_vol, zero(eltype(du_vol)))
            for i in axes(du, 1)
                for j in axes(du, 1)
                    u_i = u[i, e]
                    u_j = u[j, e]
                    du_vol[i] += 2 * Q_skew[i,j] * flux_ranocha(u_i, u_j, 1, equations)
                end
            end
        else
            du_vol .= -QTr * flux.(u_q[:, e], 1, equations)
        end
        view(du, :, e) .+= du_vol

        # determine viscosity necessary for entropy stability
        # based on local entropy error.
        psi_boundary = psi(uM[end, e], equations) - psi(uM[1, e], equations)        
        entropy_error = sum(dot.(v[:, e], du_vol)) + psi_boundary

        sigma_q = rd.Vq * sigma[:, e]
        unscaled_viscous_dissipation = zero(eltype(u[1]))
        integrand_L2_norm = zero(eltype(u[1]))
        for i in axes(sigma_q, 1)
            # v_i = cons2entropy(u_q[i, e], equations)
            v_i = cons2entropy(u_avg[i, e], equations)
            K = dudv(v_i, equations)
            unscaled_viscous_dissipation += 
                md.wJq[i, e] * (sigma_q[i]' * K * sigma_q[i])
            integrand[i] = sigma_q[i]' * K * sigma_q[i]
            integrand_L2_norm += md.wJq[i, e] * integrand[i]^2
        end
        
        # add viscosity if entropy error is positive (the wrong sign)        

        # piecewise constant version
        a = -min(0, entropy_error)
        b = unscaled_viscous_dissipation
        if params.AV_type == :elementwise
            view(viscosity, :, e) .= regularized_ratio(a, b)
        elseif params.AV_type == :subcell
            view(viscosity, :, e) .= 
                regularized_ratio.(-min(0, entropy_error) * integrand, integrand_L2_norm)
        else                 
            view(viscosity, :, e) .= zero(eltype(viscosity))
        end
    end

    # add viscous terms, scale derivs pointwise by the viscosity
    # sigma .= rd.Pq * (viscosity .* (dudv.(cons2entropy.(u_q, equations), equations)) .* (rd.Vq * sigma))
    sigma .= rd.Pq * (viscosity .* (dudv.(cons2entropy.(u_avg, equations), equations)) .* (rd.Vq * sigma))

    # compute divergence
    Q_skew = 0.5 * (rd.M * rd.Dr - rd.Dr' * rd.M)
    sigmaM = sigma[rd.Fmask, :]
    sigmaP = sigmaM[md.mapP]
    du .-= md.rxJ .* (Q_skew * sigma) + E' * (@. 0.5 * sigmaP * md.nx)

    du .= inv(rd.M) * (-du ./ md.J)
    
end

function initial_condition_density_wave(x, equations::CompressibleEulerEquations1D; 
                                        amplitude = .98)
    rho = 1.0 + amplitude * sin(2 * pi * x)
    u = .1
    p = 10

    return prim2cons(SVector(rho, u, p), equations)
end

function run_density_wave(; N=5, cells_per_dimension=4, approximation_type = SBP(),
                            amplitude = 0.5, use_EC_volume_flux = false, 
                            AV_type = :none, use_entropy_projection = true)

    rd = RefElemData(Line(), approximation_type, N)
    md = MeshData(uniform_mesh(Line(), cells_per_dimension), rd; 
                is_periodic=true) 
    u = rd.Pq * initial_condition_density_wave.(rd.Vq * md.x, equations; amplitude)

    params = (; rd, md, equations, interface_flux = flux_lax_friedrichs, 
                sigma=similar(u), viscosity=ones(size(md.xq)),
                use_EC_volume_flux, AV_type, use_entropy_projection, amplitude)

    tspan = (0.0, 25.0) # density wave
    ode = ODEProblem(rhs!, u, tspan, params)
    sol = solve(ode, SSPRK43(), 
                dt = 1e-8,
                abstol=1e-8, reltol=1e-6,
                saveat=LinRange(tspan..., 250), 
                callback=AliveCallback(alive_interval=500))

    return sol            
end


N = 7
cells_per_dimension = 4
approximation_type = SBP()
# approximation_type = Polynomial()
amplitude = 0.98
use_entropy_projection = true
sol_ESAV = run_density_wave(; N, cells_per_dimension, 
                            approximation_type, use_entropy_projection,
                            amplitude, use_EC_volume_flux = false, 
                            AV_type = :elementwise)
sol_ESAV_subcell = run_density_wave(; N, cells_per_dimension, 
                                      approximation_type, use_entropy_projection,
                                      amplitude, use_EC_volume_flux = false, 
                                      AV_type = :subcell)
sol_EC = run_density_wave(; N, cells_per_dimension, 
                            approximation_type, use_entropy_projection,
                            amplitude, use_EC_volume_flux = true, 
                            AV_type = :none)          
sol_DG = run_density_wave(; N, cells_per_dimension, 
                            approximation_type, use_entropy_projection =false,
                            amplitude, use_EC_volume_flux = false, 
                            AV_type = :none)                                      


function compute_L2_error(sol)
    (; rd, md, amplitude) = sol.prob.p
    L2_error = Float64[]
    for (i, u) in enumerate(sol.u)
        ptwise_error = rd.Vq * u - 
            initial_condition_density_wave.(md.xq .- 0.1 * sol.t[i], equations; amplitude)
        push!(L2_error, sqrt(sum(md.wJq .* norm.(ptwise_error).^2)))
    end
    return L2_error
end

plot(sol_DG.t[2:end], compute_L2_error(sol_DG)[2:end], linewidth=2, label="DG")
plot!(sol_EC.t[2:end], compute_L2_error(sol_EC)[2:end], linewidth=2, linestyle=:dash, label="EC")
plot!(sol_ESAV.t[2:end], compute_L2_error(sol_ESAV)[2:end], linewidth=2, linestyle=:dot, label="AV")
# plot!(sol_ESAV_subcell.t[2:end], compute_L2_error(sol_ESAV_subcell)[2:end], linewidth=2, linestyle=:dashdot, label="Subcell AV")
xlabel!("Time", xguidefontsize=14)
plot!(yaxis=:log, legend=:bottomleft, dpi=600, 
      xtickfontsize=14, ytickfontsize=14, legendfontsize=9)
if amplitude ≈ 0.5
    png("Density_wave_over_time_A0p5.png")
else
    png("Density_wave_over_time_A0p98.png")
end





function compute_entropy(sol)
    (; rd, md) = sol.prob.p
    entropy_over_time = [sum(md.wJq .* entropy.(rd.Vq * u, equations)) for u in sol.u]
    return entropy_over_time .- entropy_over_time[1]
end

plot(sol_DG.t, compute_entropy(sol_DG), linewidth=2, label="DG")
plot!(sol_EC.t, compute_entropy(sol_EC), linewidth=2, linestyle=:dash, label="EC")
plot!(sol_ESAV.t, compute_entropy(sol_ESAV), linewidth=2, linestyle=:dot, label="AV")
# plot!(sol_ESAV_subcell.t, compute_entropy(sol_ESAV_subcell), linewidth=2, linestyle=:dashdot, label="Subcell AV")
xlabel!("Time", xguidefontsize=14)
# ylabel!("Change in entropy", yguidefontsize=14)
# plot!(ylims=(-6e-5, 2e-5)) # for A = 0.5
plot!(xtickfontsize=12, ytickfontsize=12, legendfontsize=9, 
     legend=:bottomleft, dpi=600)
# plot!(ylims=(-7e-5, 2.5e-5))
if amplitude ≈ 0.5
    png("Density_wave_entropy_over_time_A0p5.png")
else
    png("Density_wave_entropy_over_time_A0p98.png")
end

