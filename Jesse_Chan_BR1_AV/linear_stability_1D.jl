using StartUpDG
using Plots, LinearAlgebra
using Trixi
using StaticArrays
using OrdinaryDiffEq

# loop over discretization parameters
for (N, K) in [(3, 8), (3, 16), (3, 32), (7, 4), (7, 8), (7, 16)]

rd = RefElemData(Line(), SBP(), N)
(VX,), EToV = uniform_mesh(rd.element_type, K)
md = MeshData(((VX,), EToV), rd; is_periodic=true) 

equations = CompressibleEulerEquations1D(1.4)

psi(u, ::CompressibleEulerEquations1D) = u[2]
psi(u, normal, ::CompressibleEulerEquations1D) = u[2] * normal

regularized_ratio(a, b; tol=1e-14) = a * b / (b^2 + tol)

using Trixi.ForwardDiff
dudv(v, equations) = ForwardDiff.jacobian(v -> entropy2cons(v, equations), v)

function calc_dg_derivative!(dudx, u, params)
    (; rd, md) = params
    uf = rd.Vf * u
    uP = uf[md.mapP]
    interface_flux = @. 0.5 * (uP - uf) * md.nx #- 0.5 * (uP - uf)
    dudx .= md.rxJ .* (rd.Dr * u) + rd.LIFT * interface_flux
    dudx ./= md.J
    return dudx
end

function calc_dg_derivative(u, params)
    dudx = similar(u)
    return calc_dg_derivative!(dudx, u, params)
end

function rhs!(du, u, params, t)

    (; rd, md, equations, interface_flux) = params

    # compute entropy variables
    viscosity = ones(eltype(u[1]), size(md.xq))
    v = rd.Pq * cons2entropy.(rd.Vq * u, equations)

    # calculate derivatives for viscous terms
    sigma = calc_dg_derivative(v, params)

    u_q = rd.Vq * u
    u_avg = repeat(0.5 * sum(rd.M, dims=1) * u, length(rd.rq), 1)

    # solution "u" to evaluate the AV at
    # u_av = u_avg
    u_av = u_q
    # u_av = entropy2cons.(rd.Vq * v, equations)

    uM = entropy2cons.(rd.Vf * v, equations)
    uP = uM[md.mapP]

    if md.is_periodic[1] == false
        uP[1] = initial_condition(md.x[1], equations)
        uP[end] = initial_condition(md.x[end], equations)
    end

    # surface terms
    E = rd.Vf   
    du .= E' * (@. interface_flux(uM, uP, SVector(md.nx), equations))

    # volume terms
    QTr = rd.Dr' * rd.M * rd.Pq
    Q_skew = 0.5 * (rd.M * rd.Dr - rd.Dr' * rd.M)

    du_vol = similar(du[:, 1])
    integrand = zeros(eltype(u[1]), size(u_q, 1))
    for e in axes(du, 2)
        if rd.approximation_type isa SBP
            fill!(du_vol, zero(eltype(du_vol)))
            for i in axes(u_q, 1), j in axes(u_q, 1)
                if i > j
                    u_i = u_q[i, e]
                    u_j = u_q[j, e]
                    Fij = 2 * Q_skew[i, j] * params.volume_flux(u_i, u_j, 1, equations)
                    du_vol[i] += Fij
                    du_vol[j] -= Fij
                end
            end
        else # if modal
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
            v_i = cons2entropy(u_av[i, e], equations)
            K = dudv(v_i, equations)
            unscaled_viscous_dissipation += 
                md.wJq[i, e] * (sigma_q[i]' * K * sigma_q[i])
            integrand[i] = sigma_q[i]' * K * sigma_q[i]
            integrand_L2_norm += md.wJq[i, e] * integrand[i]^2
        end
        
        # piecewise constant viscosity
        # smoothmin(a, b; alpha=100) = (a * exp(-alpha * a) + b * exp(-alpha * b)) / (exp(-alpha * a) + exp(-alpha * b))
        # a = -smoothmin(0, entropy_error)
        a = -min(0, entropy_error)
        b = unscaled_viscous_dissipation
        view(viscosity, :, e) .= regularized_ratio(a, b)

        # # # subcell viscosity
        # view(viscosity, :, e) .= 
        #     regularized_ratio.(-min(0, entropy_error) * integrand, integrand_L2_norm)

    end

    # add viscous terms, scale derivs pointwise by the viscosity
    sigma .= rd.Pq * (viscosity .* (dudv.(cons2entropy.(u_av, equations), equations)) .* (rd.Vq * sigma))

    # compute divergence
    sigmaM = sigma[rd.Fmask, :]
    sigmaP = sigmaM[md.mapP]
    if params.use_AV == true
        du .-= md.rxJ .* (Q_skew * sigma) + E' * (@. 0.5 * sigmaP * md.nx) #+ 0.5 * (sigmaP - sigmaM))
    end

    du .= inv(rd.M) * (-du ./ md.J)
    return du
end

# rho = ..., rho * u = ..., rho * e = ...
function initial_condition(x, equations::CompressibleEulerEquations1D; amplitude=.98)

    rho = 1.0 + amplitude * sin(2 * pi * x)
    u = .10
    p = 10.0

    return prim2cons(SVector(rho, u, p), equations)
end

# initialize using L2 projection
rd2 = RefElemData(Line(), rd.N) #; quad_rule_vol=gauss_lobatto_quad(0,0,rd.N+4))
u = rd2.Pq * initial_condition.(rd2.Vq * md.x, equations; amplitude=0.5)
# u = initial_condition.(md.x, equations)

interface_flux = flux_lax_friedrichs
# interface_flux = flux_central

volume_flux = flux_central
# volume_flux = flux_ranocha

use_AV = true
use_AV = false

params = (; rd, md, equations, interface_flux, volume_flux, use_AV)

tspan = (0, 2.0)
ode = ODEProblem(rhs!, u, tspan, params)
sol = solve(ode, SSPRK43(), 
            dt = 1e-8, abstol=1e-8, reltol=1e-6,
            #callback=AliveCallback(alive_interval=1000))
            )
u = sol.u[end]

using ForwardDiff
rhs_fd = let params=params
    function rhs_fd(u_in::AbstractArray{T}) where {T}
        u = reinterpret(reshape, SVector{nvariables(equations), T}, 
                        reshape(u_in, :, size(md.x)...))
        du = fill!(similar(u), zero(eltype(u)))
        rhs!(du, u, params, 0.0)
        return reinterpret(reshape, T, du)
    end
end

u_in = reinterpret(reshape, Float64, u)
A = ForwardDiff.jacobian(rhs_fd, u_in)

# # compute Jacobian using finite differences instead
# duP, duM = ntuple(_ -> similar(u), 2)
# delta_u = zeros(3, size(md.x)...)
# du = similar(u)
# A = zeros(length(delta_u), length(delta_u))
# tol = 1e-7
# for i in eachindex(delta_u)
#     fill!(du, zero(eltype(du)))
#     delta_u[i] = tol
#     rhs!(duP, u .+ reinterpret(reshape, SVector{3, eltype(delta_u)}, delta_u), params, 0.0)
#     rhs!(duM, u .- reinterpret(reshape, SVector{3, eltype(delta_u)}, delta_u), params, 0.0)
#     A[:,i] = vec(reinterpret(reshape, eltype(delta_u), (duP - duM)) ./ (2 * tol))
#     delta_u[i] = 0
# end

@show extrema(real.(eigvals(A)))

# scatter(eigvals(A), leg=false, ms=6)#, label="AV", marker=:star5)
# plot!(xlims=(-5, 5), ylims=(-5, 5))
# plot!(xlims=(-.25, .25), ylims=(-5, 5))
# plot!(dpi=500, tickfontsize=14)
# xlabel!(""); ylabel!("")

end # loop over (N, K)