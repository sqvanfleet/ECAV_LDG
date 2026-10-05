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



N = 3
rd = RefElemData(Line(), N)
(VX,), EToV = uniform_mesh(rd.element_type, 100)

enforce_additional_entropy_inequality = false
use_EC_volume_flux = false

# @. VX = 0.5 * (1 + VX) # modified Sod
# init_condition = initial_condition_modified_sod
# tspan = (0.0, 0.2) # modified Sod

@. VX = 5 * VX # for Shu Osher
init_condition = initial_condition_shu_osher
tspan = (0.0, 1.8) # Shu-Osher

md = MeshData(((VX,), EToV), rd)

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
# interface_flux = flux_ranocha
# interface_flux = flux_central

params_elementwise = (; rd, md, equations, interface_flux, init_condition,
            sigma=similar(u), viscosity=ones(size(md.xq)),
            enforce_additional_entropy_inequality, 
            use_EC_volume_flux, use_entropy_projection = true,
            AV_type = :elementwise)      

ode = ODEProblem(rhs!, u, tspan, params_elementwise)
sol_elementwise = solve(ode, SSPRK43(), 
                        dt = 1e-8,
                        abstol=1e-8, reltol=1e-6,
                        saveat=LinRange(tspan..., 100), 
                        callback=AliveCallback(alive_interval=100))

params_subcell = (; rd, md, equations, interface_flux, init_condition,
            sigma=similar(u), viscosity=ones(size(md.xq)),
            enforce_additional_entropy_inequality, 
            use_EC_volume_flux, use_entropy_projection = true, 
            AV_type = :subcell)    

ode = ODEProblem(rhs!, u, tspan, params_subcell)
sol_subcell = solve(ode, SSPRK43(), 
                    dt = 1e-8,
                    abstol=1e-8, reltol=1e-6,
                    # abstol=1e-9, reltol=1e-7,
                    saveat=LinRange(tspan..., 100), 
                    callback=AliveCallback(alive_interval=100))            

pad_nans(u) = vec([u; fill(NaN, 1, size(u, 2))])

# using LaTeXStrings
# plot(pad_nans(rd.Vp * md.x), pad_nans(rd.Vp * getindex.(u, 1)), linewidth=2,
#       label=L"N=%$N, K=%$(md.num_elements)")
# plot!(legend=:topright, dpi=500, xformatter=:none, yformatter=:none, legendfontsize=14)

using MAT
vars = matread("weno5_shuosher.mat")   
plot(vec(vars["x"]), vec(vars["rho"]), label="WENO", ylims=(0.5, 5.5), linewidth=1.5)

Vp = vandermonde(Line(), rd.N, LinRange(-1, 1, 100)) / rd.VDM

u = sol_elementwise.u[end]
plot!(pad_nans(Vp * md.x), pad_nans(Vp * getindex.(u, 1)), linewidth=2, 
      label=L"$N=%$N, K=%$(md.num_elements)$, elementwise AV")

u = sol_subcell.u[end]
plot!(pad_nans(Vp * md.x), pad_nans(Vp * getindex.(u, 1)), linewidth=2,
      label=L"$N=%$N, K=%$(md.num_elements)$, subcell AV")


plot!(legend=:bottomleft, xformatter=:none, yformatter=:none, dpi=500)
plot!(dpi=500, xformatter=:none, yformatter=:none, legendfontsize=14)
#
plot!(xlims=(0.5, 2.5), ylims=(2.25, 4.8)) # zoom 
# png("shu_osher_elementwise_vs_subcell.png")

