using StartUpDG
using Plots, LinearAlgebra
using Trixi
using StaticArrays
using OrdinaryDiffEq

function initial_condition_density_wave(x, equations::CompressibleEulerEquations1D; 
                                        amplitude = .98)
    rho = 1.0 + amplitude * sin(2 * pi * x)
    u = .1
    p = 10

    return prim2cons(SVector(rho, u, p), equations)
end

# modified Sod
function initial_condition_modified_sod(x, equations::CompressibleEulerEquations1D)
    if x < 0.3
        return prim2cons(SVector(1, .75, 1), equations)
    else
        return prim2cons(SVector(0.125, 0.0, .1), equations)
        # return prim2cons(SVector(0.0125, 0.0, .01), equations)
    end
end

# Shu-Osher
function initial_condition_shu_osher(x, equations::CompressibleEulerEquations1D)
    if x < -4
        rho = 3.857143        
        u = 2.629369
        p = 10.3333
    else
        rho = 1 + .2 * sin(5 * x)
        u = 0.0
        p = 1
    end
    return prim2cons(SVector(rho, u, p), equations)        
end



N = 7
rd = RefElemData(Line(), N)
(VX,), EToV = uniform_mesh(rd.element_type, 50)

enforce_additional_entropy_inequality = false
use_EC_volume_flux = false
AV_type = :elementwise

@. VX = 0.5 * (1 + VX) # modified Sod
init_condition = initial_condition_modified_sod
tspan = (0.0, 0.2) # modified Sod

# @. VX = 5 * VX # for Shu Osher
# init_condition = initial_condition_shu_osher
# tspan = (0.0, 1.8) # Shu-Osher

md = MeshData(((VX,), EToV), rd)

# # for density and blast wave
# init_condition = initial_condition_density_wave
# tspan = (0.0, 10.0) # density wave
# md = MeshData(((VX,), EToV), rd; is_periodic=true) 

equations = CompressibleEulerEquations1D(1.4)

psi(u, ::CompressibleEulerEquations1D) = u[2]
psi(u, normal, ::CompressibleEulerEquations1D) = u[2] * normal

# adapted from Atum.jl
function EC_matrix_dissipation(u_ll, u_rr, n⃗::SVector, 
                               equations::CompressibleEulerEquations1D)
    
    γ = equations.gamma
    ecflux = flux_ranocha(u_ll, u_rr, n⃗, equations)

    # in 1D, n is just a scalar
    n⃗ = n⃗[1]

    ρ⁻, ρu⃗⁻, ρe⁻ = u_ll
    u⃗⁻ = ρu⃗⁻ / ρ⁻
    e⁻ = ρe⁻ / ρ⁻
    p⁻ = Trixi.pressure(u_ll, equations)
    b⁻ = ρ⁻ / 2p⁻

    ρ⁺, ρu⃗⁺, ρe⁺ = u_rr
    u⃗⁺ = ρu⃗⁺ / ρ⁺
    e⁺ = ρe⁺ / ρ⁺
    p⁺ = Trixi.pressure(u_rr, equations)
    b⁺ = ρ⁺ / 2p⁺

    logavg = Trixi.ln_mean
    avg(a, b) = 0.5 * (a + b)
    ρ_log = logavg(ρ⁻, ρ⁺)
    b_log = logavg(b⁻, b⁺)
    u⃗_avg = avg(u⃗⁻, u⃗⁺)
    p_avg = avg(ρ⁻, ρ⁺) / 2avg(b⁻, b⁺)
    u²_bar = 2 * norm(u⃗_avg) - avg(norm(u⃗⁻), norm(u⃗⁺))
    h_bar = γ / (2 * b_log * (γ - 1)) + u²_bar / 2
    c_bar = sqrt(γ * p_avg / ρ_log)

    u⃗mc = u⃗_avg - c_bar * n⃗
    u⃗pc = u⃗_avg + c_bar * n⃗
    u_avgᵀn = u⃗_avg * n⃗

    v⁻ = cons2entropy(u_ll, equations)
    v⁺ = cons2entropy(u_rr, equations)
    #v⁻ = Atum.entropyvariables(law, q⁻, aux⁻)
    #v⁺ = Atum.entropyvariables(law, q⁺, aux⁺)
    Δv = v⁺ - v⁻

    λ1 = abs(u_avgᵀn - c_bar) * ρ_log / 2γ
    λ2 = abs(u_avgᵀn) * ρ_log * (γ - 1) / γ
    λ3 = abs(u_avgᵀn + c_bar) * ρ_log / 2γ
    λ4 = abs(u_avgᵀn) * p_avg

    Δv_ρ, Δv_ρu⃗, Δv_ρe = Δv
    u⃗ₜ = u⃗_avg - u_avgᵀn * n⃗

    w1 = λ1 * (Δv_ρ + u⃗mc' * Δv_ρu⃗ + (h_bar - c_bar * u_avgᵀn) * Δv_ρe)
    w2 = λ2 * (Δv_ρ + u⃗_avg' * Δv_ρu⃗ + u²_bar / 2 * Δv_ρe)
    w3 = λ3 * (Δv_ρ + u⃗pc' * Δv_ρu⃗ + (h_bar + c_bar * u_avgᵀn) * Δv_ρe)

    Dρ = w1 + w2 + w3

    Dρu⃗ = (w1 * u⃗mc +
           w2 * u⃗_avg +
           w3 * u⃗pc +
           λ4 * (Δv_ρu⃗ - n⃗' * (Δv_ρu⃗) * n⃗ + Δv_ρe * u⃗ₜ))

    Dρe = (w1 * (h_bar - c_bar * u_avgᵀn) +
           w2 * u²_bar / 2 +
           w3 * (h_bar + c_bar * u_avgᵀn) +
           λ4 * (u⃗ₜ' * Δv_ρu⃗ + Δv_ρe * (u⃗_avg' * u⃗_avg - u_avgᵀn ^ 2)))

    return ecflux - SVector(Dρ, Dρu⃗..., Dρe) / 2
end

# use the 2D HLLC implementation since the 1D Trixi version doesn't
# account for n::Integer < 0.
import Trixi: flux_hllc
function flux_hllc(u_ll, u_rr, n::SVector{1}, equations::CompressibleEulerEquations1D)
    f = flux_hllc(SVector(u_ll[1], u_ll[2], 0, u_ll[3]), 
                  SVector(u_rr[1], u_rr[2], 0, u_rr[3]), 
                  SVector(n[1], 0.0), 
                  CompressibleEulerEquations2D(equations.gamma))    
    return SVector(f[1], f[2], f[4])
end

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

    if md.is_periodic[1] == false
        (; init_condition) = params
        uP[1] = init_condition(md.x[1], equations)
        uP[end] = init_condition(md.x[end], equations)
    end

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

        # enforce the max between the ESAV entropy inequality and the 
        # entropy inequality where psi(u) is evaluated with u_h instead 
        # of the entropy projection ũ = u(P * v(u_h))
        if params.enforce_additional_entropy_inequality == true
            # entropy error where psi is evaluated using u_h instead 
            # of the entropy projection. 
            entropy_error_2 = 
                sum(dot.(v[:, e], du_vol)) + 
                psi(u[end, e], equations) - psi(u[1, e], equations)

            a = max(a, -min(0, entropy_error_2))
        end

        b = unscaled_viscous_dissipation
        
        if params.AV_type == :elementwise
            view(viscosity, :, e) .= regularized_ratio(a, b)
        elseif params.AV_type == :subcell
            view(viscosity, :, e) .= 
                regularized_ratio.(-min(0, entropy_error) * integrand, integrand_L2_norm)
        else                 
            view(viscosity, :, e) .= zero(eltype(viscosity))
        end
        
        # view(viscosity, :, e) .= regularized_ratio(a, b)

        # # subcell version
        # view(viscosity, :, e) .= 
        #     regularized_ratio.(-min(0, entropy_error) * integrand, integrand_L2_norm)

    end

    # # # enforce C0 continuity between elements    
    # visc_polynomial = rd.Pq * viscosity
    # viscosity_M = rd.Vf * visc_polynomial
    # viscosity_P = viscosity_M[md.mapP]
    # viscosity .= rd.Vq * rd.V1 * max.(viscosity_M, viscosity_P)

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

u = rd.Pq * init_condition.(rd.Vq * md.x, equations)

interface_flux = flux_lax_friedrichs
# interface_flux = flux_hllc
# interface_flux = EC_matrix_dissipation

params = (; rd, md, equations, interface_flux, init_condition,
            sigma=similar(u), viscosity=ones(size(md.xq)),
            enforce_additional_entropy_inequality, 
            use_EC_volume_flux, AV_type, 
            use_entropy_projection = true)


ode = ODEProblem(rhs!, u, tspan, params)
sol = solve(ode, SSPRK43(), 
            dt = 1e-8,
            abstol=1e-8, reltol=1e-6,
            saveat=LinRange(tspan..., 100), 
            callback=AliveCallback(alive_interval=100))

u = sol.u[end]
pad_nans(u) = vec([u; fill(NaN, 1, size(u, 2))])

using LaTeXStrings
plot(pad_nans(rd.Vp * md.x), pad_nans(rd.Vp * getindex.(u, 1)), linewidth=2,
      label=L"N=%$N, K=%$(md.num_elements)")
plot!(legend=:topright, dpi=500, xformatter=:none, yformatter=:none, legendfontsize=14)

# using MAT
# vars = matread("weno5_shuosher.mat")   
# plot(vec(vars["x"]), vec(vars["rho"]), label="WENO", ylims=(0.5, 5.5), linewidth=2)
# plot!(pad_nans(rd.Vp * md.x), pad_nans(rd.Vp * getindex.(u, 1)), linewidth=2,
#       label=L"N=%$N, K=%$(md.num_elements)")
# plot!(legend=:bottomleft, xformatter=:none, yformatter=:none, dpi=500)
# plot!(dpi=500, xformatter=:none, yformatter=:none, legendfontsize=14)
# plot!(xlims=(-4, 3), ylims=(2.5, 5)) # zoom 


# # using LaTeXStrings
# # if params.enforce_additional_entropy_inequality==true
# #     plot!(pad_nans(rd.Vp * md.x), pad_nans(rd.Vp * getindex.(u, 1)),       
# #         label="Two entropy inequalities", linewidth=2)
# # else
# #     plot!(pad_nans(rd.Vp * md.x), pad_nans(rd.Vp * getindex.(u, 1)),       
# #         label="One entropy inequality", linewidth=2)
# # end
# plot!(pad_nans(rd.Vp * md.x), pad_nans(rd.Vp * getindex.(u, 1)),       
#       label=L"N=%$N, K=%$(md.num_elements)")
# plot!(legend=:topleft, xformatter=:none, yformatter=:none, dpi=500)
# plot!(dpi=500, xformatter=:none, yformatter=:none)      