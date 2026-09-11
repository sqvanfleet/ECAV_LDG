using OrdinaryDiffEq
using StartUpDG
using Plots, LinearAlgebra
using Trixi
using StaticArrays
using Polyester, Octavian
using StructArrays

function dudv_cons_vars(u, equations::CompressibleEulerEquations2D)     

    _, rho_v1, rho_v2, E = u
    rho, v1, v2, p = cons2prim(u, equations)
    
    A11, A12, A13, A14 = rho, rho_v1, rho_v2, E
    A22 = rho_v1 * v1 + p
    A23 = rho_v1 * v2
    A24 = v1 * (E + p)
    A33 = rho_v2 * v2 + p
    A34 = v2 * (E + p)
    
    # internal energy
    a2 = equations.gamma * p / rho
    H = a2 * equations.inv_gamma_minus_one + 0.5 * (v1^2 + v2^2)
    A44 = rho * H * H - a2 * p * equations.inv_gamma_minus_one

    return SMatrix{4, 4, eltype(u), 16}(A11, A12, A13, A14, 
                                        A12, A22, A23, A24,
                                        A13, A23, A33, A34,
                                        A14, A24, A34, A44)
end

psi(u, normal, ::CompressibleEulerEquations2D) = u[2] * normal[1] + u[3] * normal[2]

function mymul!(du::StructArray, A, u, alpha=1, beta=0)

    # return LinearAlgebra.mul!(du, A, u, alpha, beta)

    # mul_by! = (du, u) -> LinearAlgebra.mul!(du, A, u, alpha, beta)

    # du .= beta * du + alpha * A * u
    # return du

    mul_by! = (du, u) -> Octavian.matmul!(du, A, u, alpha, beta)
    StructArrays.foreachfield(mul_by!, du, u)
    return du
end

regularized_ratio(a, b) = a * b / (b^2 + 1e-14)

function calc_dg_gradients!(gradients, u, params)
    (; rd, md) = params
    dudx = fill!(gradients[1], zero(eltype(gradients[1])))
    dudy = fill!(gradients[2], zero(eltype(gradients[2])))

    (; uf, interface_flux_x, interface_flux_y) = params
    mymul!(uf, rd.Vf, u)
    @batch for i in eachindex(interface_flux_x, interface_flux_y)
        interface_flux_x[i] = 0.5 * (uf[md.mapP[i]] - uf[i]) * md.nxJ[i]
        interface_flux_y[i] = 0.5 * (uf[md.mapP[i]] - uf[i]) * md.nyJ[i]
    end

    (; dudr, duds) = params
    mymul!(dudr, rd.Dr, u)
    mymul!(duds, rd.Ds, u)
    mymul!(dudx, rd.LIFT, interface_flux_x)
    mymul!(dudy, rd.LIFT, interface_flux_y)
    @batch for i in eachindex(dudx)
        invJ = inv(md.J[i])
        dudx[i] += (md.rxJ[i] * dudr[i] + md.sxJ[i] * duds[i]) * invJ
        dudy[i] += (md.ryJ[i] * dudr[i] + md.syJ[i] * duds[i]) * invJ
    end
end

function rhs!(du, u, params, t)

    (; rd, md, equations, interface_flux) = params
    (; invMQTr, invMQTs, invMQr_skew, invMQs_skew) = params
    (; K_visc) = params

    # interpolate to quad points, project entropy variables
    uq = rd.Vq * u
    vq = similar(uq)
    @batch for i in eachindex(vq)
        vq[i] = cons2entropy(uq[i], equations)
    end
    v = rd.Pq * vq

    # calculate derivatives for viscous terms
    calc_dg_gradients!(params.gradients, v, params)

    # interpolate
    vf = rd.Vf * v   
    uM = vf # alias
    @batch for i in eachindex(uM) 
        uM[i] = entropy2cons(vf[i], equations)
    end
    uP = uM[md.mapP]
    if all(md.is_periodic) == false
        @. uP[md.mapB] = initial_condition(md.xf[md.mapB], md.yf[md.mapB], t, equations)
    end

    # compute volume contributions
    fill!(du, zero(eltype(du)))
    (; rxyJ, sxyJ, flux_uq) = params
    @batch for i in eachindex(uq)
        flux_uq[i] = flux(uq[i], rxyJ[i], equations)
    end
    mymul!(du, invMQTr, flux_uq)
    @batch for i in eachindex(uq)
        flux_uq[i] = flux(uq[i], sxyJ[i], equations)
    end
    mymul!(du, invMQTs, flux_uq, 1, 1)

    # compute minimal artificial viscosity coefficients
    psi_boundary = sum(Diagonal(rd.wf) * (@. psi(uM, SVector.(md.nxJ, md.nyJ), equations)), dims=1)
    entropy_errors = sum(dot.(v, rd.M * du), dims=1) + psi_boundary

    (; sigma_x, sigma_y) = params
    mul!(sigma_x, rd.Vq, params.gradients[1])
    mul!(sigma_y, rd.Vq, params.gradients[2])

    u_avg = params.eTrM * u
    @batch for e in axes(K_visc, 2)        
        for i in axes(K_visc, 1)
            K_visc[i, e] = dudv_cons_vars(u_avg[e], equations)
        end
    end

    if params.av_type != :subcell
        # element-local K_visc                   
        unscaled_viscous_dissipation = zeros(md.num_elements)
        @batch for e in eachindex(unscaled_viscous_dissipation)
            val = zero(eltype(unscaled_viscous_dissipation))
            for i in axes(md.wJq, 1)
                val += md.wJq[i, e] * (dot(sigma_x[i, e], K_visc[i, e], sigma_x[i, e]) + 
                                    dot(sigma_y[i, e], K_visc[i, e], sigma_y[i, e]))
            end
            unscaled_viscous_dissipation[e] = val
        end
        (; viscosity) = params
        viscosity .= vec(map((a,b) -> max(0.0, -min(0, a) * b / (1e-12 + b * b)), 
                        entropy_errors, unscaled_viscous_dissipation))
        @batch for e in axes(K_visc, 2)
            for i in axes(K_visc, 1)
                K_visc[i, e] *= viscosity[e]
            end
        end
    elseif params.av_type == :subcell
        # subcell varying viscosity
        (; viscosity, integrand) = params
        fill!(viscosity, zero(eltype(viscosity)))
        @batch for e in axes(K_visc, 2)
            integrand_L2_norm = 0.0
            for i in eachindex(integrand)
                integrand[i] = (dot(sigma_x[i, e], K_visc[i, e], sigma_x[i, e]) + 
                                dot(sigma_y[i, e], K_visc[i, e], sigma_y[i, e]))
                integrand_L2_norm += md.wJq[i, e] * integrand[i]^2
            end
            # rescale K_visc
            for i in axes(K_visc, 1)
                subcell_viscosity = 
                    regularized_ratio(-min(0, entropy_errors[e]) * integrand[i], 
                                    integrand_L2_norm)
                params.subcell_viscosity[i, e] = subcell_viscosity
                K_visc[i, e] *= subcell_viscosity
            end
        end

    else # if not specified, set to zero
        K_visc .= zero(eltype(K_visc))
    end # if 


    # now add transport surface terms
    (; convective_interface_flux) = params
    @batch for i in eachindex(convective_interface_flux)
        convective_interface_flux[i] = 
            interface_flux(uM[i], uP[i], SVector(md.nx[i], md.ny[i]), equations) * md.Jf[i]
    end
    mymul!(du, rd.LIFT, convective_interface_flux, 1, 1)

    # compute scaled gradients and their divergence
    @batch for i in eachindex(sigma_x)
        sigma_x[i] = K_visc[i] * sigma_x[i]
        sigma_y[i] = K_visc[i] * sigma_y[i]
    end
    (; sigma_x_projected, sigma_y_projected) = params
    mymul!(sigma_x_projected, rd.Pq, sigma_x)
    mymul!(sigma_y_projected, rd.Pq, sigma_y)

    (; sigma_r, sigma_s, sigma_x_f, sigma_y_f, sigmaM) = params
    mymul!(sigma_x_f, rd.Vf, sigma_x_projected)
    mymul!(sigma_y_f, rd.Vf, sigma_y_projected)
    @. sigmaM = md.nxJ * sigma_x_f + md.nyJ * sigma_y_f
    @batch for e in axes(sigma_x_projected, 2)
        for i in axes(sigma_x_projected, 1)
            sigma_r[i, e] = md.rxJ[1, e] * sigma_x_projected[i, e] + md.ryJ[1, e] * sigma_y_projected[i, e]
            sigma_s[i, e] = md.sxJ[1, e] * sigma_x_projected[i, e] + md.syJ[1, e] * sigma_y_projected[i, e]
        end
    end

    # sigmaM[md.mapP] is negated since sigmaM is scaled by the normal
    (; sigma_interface_flux) = params
    @. sigma_interface_flux = -0.5 * sigmaM[md.mapP]
    mymul!(du, invMQr_skew, sigma_r, -1, 1)
    mymul!(du, invMQs_skew, sigma_s, -1, 1)
    mymul!(du, rd.LIFT, sigma_interface_flux, -1, 1)
    
    @batch for i in eachindex(du)
        du[i] /= -md.J[i]
    end
    
end

# # rho = ..., rho * u = ..., rho * e = ...
function initial_condition_density_wave(x, y, t, equations::CompressibleEulerEquations2D)
    u = .10 
    v = .20
    rho = 1.0 + .5 * sin(2.0 * pi * ((x - u * t) + (y - v * t)))
    # rho = 1.0 + .5 * sin(pi * (x - u * t)) * sin(pi * (y - v * t))

    p = 10.0

    return prim2cons(SVector(rho, u, v, p), equations)
end

function compute_L2_error(N, cells_per_dimension)
    rd = RefElemData(Tri(), N)
    md = MeshData(uniform_mesh(Tri(), cells_per_dimension), rd; 
                is_periodic = true)

    equations = CompressibleEulerEquations2D(1.4)
    init_condition = initial_condition_density_wave

    u = StructArray{SVector{nvariables(equations), Float64}}(ntuple(_ -> similar(md.x), 4))
    u .= rd.Pq * init_condition.(rd.Vq * md.x, rd.Vq * md.y, 0.0, equations)

    u_plot = getindex.(parent(u),1)

    # p = scatter(vec(rd.Vp * md.x), vec(rd.Vp * md.y), zcolor = vec(rd.Vp * u_plot), 
    #         msw = 0, ms=1, title="Initial condition", 
    #             xlabel="x", ylabel="y", legend = false, cam = (0,90))

    # display(p)

    interface_flux = flux_lax_friedrichs
    # av_type = :subcell
    av_type = :elementwise
    # av_type = :none

    params = (; rd, md, equations, interface_flux, 
                v=similar(u), 
                gradients = (similar(u), similar(u)),
                viscosity=ones(size(u, 2)), 
                subcell_viscosity=ones(size(md.xq)), 
                invMQTr = rd.M \ (-rd.Dr' * rd.M * rd.Pq),
                invMQTs = rd.M \ (-rd.Ds' * rd.M * rd.Pq),
                invMQr_skew = rd.M \ (0.5 * (rd.M * rd.Dr - rd.Dr' * rd.M)),
                invMQs_skew = rd.M \ (0.5 * (rd.M * rd.Ds - rd.Ds' * rd.M)),
                rxyJ = SVector.(rd.Vq * md.rxJ, rd.Vq * md.ryJ),
                sxyJ = SVector.(rd.Vq * md.sxJ, rd.Vq * md.syJ),
                K_visc = dudv_cons_vars.((rd.Vq * u), equations),
                dudr = similar(u), duds = similar(u),
                uf = similar(rd.Vf * u),
                interface_flux_x = similar(rd.Vf * u),
                interface_flux_y = similar(rd.Vf * u),
                convective_interface_flux = similar(rd.Vf * u),            
                viscosity_quadrature = similar(md.xq),
                visc_local = zeros(3, md.num_elements), 
                flux_uq = similar(rd.Vq * u), 
                sigma_x = similar(rd.Vq * u), sigma_y = similar(rd.Vq * u),
                sigma_interface_flux = similar(rd.Vf * u),
                sigma_x_projected = similar(u),
                sigma_y_projected = similar(u),
                sigma_r = similar(u), sigma_s = similar(u),        
                eTrM = 0.5 * sum(rd.M, dims=1), 
                integrand = zeros(size(md.xq, 1)),
                sigma_x_f = similar(rd.Vf * u), 
                sigma_y_f = similar(rd.Vf * u), 
                sigmaM = similar(rd.Vf * u),
                av_type
            )

    tspan = (0.0, 1.7)

    ode = ODEProblem(rhs!, u, tspan, params)

    sol = solve(ode, SSPRK43(), 
                dt = 1e-8,
                abstol=1e-8, reltol=1e-6,
                saveat=LinRange(tspan..., 5), 
                callback=AliveCallback(alive_interval=100))

    function calc_L2_error(u, t)
        ptwise_error = rd.Vq * u - 
            initial_condition_density_wave.(md.xq, md.yq, t, equations)                                
        L2_error = sqrt(sum(md.wJq .* norm.(ptwise_error).^2))
        return L2_error
    end

    # return sol

    @show calc_L2_error(sol.u[1], sol.t[1])
    @show calc_L2_error(sol.u[end], sol.t[end])
    return calc_L2_error(sol.u[end], sol.t[end])
end

L2_errors = [compute_L2_error(N, cells_per_dimension) 
    for cells_per_dimension in [4, 8, 16, 32, 64], N in 2:4]

L2_rates = log.(L2_errors[2:end, :] ./ L2_errors[1:end-1,:]) / log(0.5)
error_table = hcat([hcat(L2_errors[:,i], [NaN; L2_rates[:,i]])  for i in axes(L2_errors, 2)]...)

using LatexPrint
using Printf
tabular([(@sprintf("%.3e", err)) for err in error_table])