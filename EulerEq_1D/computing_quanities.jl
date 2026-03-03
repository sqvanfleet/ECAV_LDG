using StartUpDG
using LinearAlgebra
using OrdinaryDiffEq, RecursiveArrayTools
using JLD2
using Trixi
using StaticArrays

#Compute the error at each saved time
function compute_L2_error_evolution_stationary_contact_wave(sol, equations, md, rd)
    L2_error = zeros(length(sol.t))
    local rq, wq = quad_nodes(rd.element_type, rd.N + 5)
    local Vq = vandermonde(rd.element_type, rd.N, rq) / rd.VDM
    local xq = Vq * md.x
    wJq = wq.*Vq*md.J
        for (i, t) in enumerate(sol.t)
            error = Vq * parent(sol.u[i]) - stationary_contact_wave.(xq, equations)
            L2_error[i] = sqrt(sum(wJq .* norm.(error).^2))
        end
    return L2_error
end

# Compute the ||dv||^2/||Pdv||^2
function compute_dvPdv(sol, equations,md, rd)
    dvPdv = zeros(length(sol.t))
    maxdvPdv_ind = 1
    for (i,u) in enumerate(sol.u)
        v = cons2entropy.(rd.Vq * parent(u),equations)
        v_avg = repeat(0.5 * sum(Diagonal(rd.wq) * (rd.Vq * rd.Pq * v), dims=1), rd.Np, 1)
        num = sum(md.wJq .* norm.(v - rd.Vq * v_avg).^2, dims = 1)
        den = sum(md.wJq .* norm.(rd.Vq * rd.Pq * (v - rd.Vq * v_avg)).^2, dims = 1)
        maxdvPdv_ind = argmax(num./den)
        dvPdv[i] = maximum(num./(den .+ 1e-14))
    end
    return dvPdv, maxdvPdv_ind
end

function compute_quantities(u_voa, rd, md, equations, AV_type)

    u = parent(u_voa)
    du = similar(u)

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
    # vM = rd.Vf * v # nodal
    # vP = vM[md.mapP] #nodal
    vM = cons2entropy.(uM, equations) #modal
    vP = cons2entropy.(uP, equations) #modal
    interface_flux = @. params.interface_flux(uM, uP, SVector(md.nxJ), equations)
    duvol = md.rxJ .* (invMQTr * rd.Pq * flux.(rd.Vq * u, 1, equations)) 
    du .= duvol + rd.LIFT * interface_flux
    if AV_type == :BR1
        beta = 0.0
    elseif AV_type == :LDG
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
    sigma = rd.Pq * (K .* (rd.Vq * theta))
    den = sum(md.wJq .* dot.(rd.Vq * sigma, rd.Vq * theta), dims=1)
    epsilon = @. num * den / (1e-14 + den^2)

    return epsilon

end

function EC_matrix_dissipation(u_ll, u_rr, n⃗::SVector, 
                               equations::CompressibleEulerEquations1D)
    
    γ = equations.gamma
    ecflux = flux_ranocha(u_ll, u_rr, n⃗, equations)

    # in 1D, n is just a scalar
    n⃗ = n⃗[1]

    ρ⁻, ρu⃗⁻, ρe⁻ = u_ll
    u⃗⁻ = ρu⃗⁻ / ρ⁻
    p⁻ = Trixi.pressure(u_ll, equations)
    b⁻ = ρ⁻ / 2p⁻

    ρ⁺, ρu⃗⁺, ρe⁺ = u_rr
    u⃗⁺ = ρu⃗⁺ / ρ⁺
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
