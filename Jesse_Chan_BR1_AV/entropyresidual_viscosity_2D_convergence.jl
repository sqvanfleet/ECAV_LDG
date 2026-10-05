using StartUpDG
using Plots, LinearAlgebra
using Trixi
using StaticArrays
using Polyester, Octavian
using StructArrays

N = 2

function initial_condition_test(x, y, t, equations::CompressibleEulerEquations2D)
    # smooth
    rho = 1 + .5 * sin(.1 + pi * x) * sin(.2 + pi * y)
    v1 = .5 * sin(.2 + pi * x) * sin(.1 + pi * y)
    v2 = 0
    p = rho^equations.gamma

    # smooth
    rho = 1 + .5 * sin(.1 + pi * x) * sin(.2 + 1.5 * pi * y)
    v1 = .5 * sin(.2 + 1.5 * pi * x) * sin(.1 + pi * y)
    v2 = 0
    p = rho^equations.gamma

    # # discontinuous
    # rho = 1.0 + (abs(0.3 * x + y) < 0.5) 
    # v1 = .1 * (abs(0.3 * x + y) < 0.5)
    # v2 = .2 * (abs(0.3 * x + y) < 0.5)
    # p = rho^equations.gamma

    return prim2cons(SVector(rho, v1, v2, p), equations)
end

max_elementwise_visc = Float64[]
max_entropy_errors = Float64[]
max_visc_dissipation = Float64[]
min_visc_dissipation = Float64[]
for K in [4, 8, 16, 32, 64]
    rd = RefElemData(Tri(), N; 
                     quad_rule_vol=quad_nodes(Tri(), N + 1), 
                     quad_rule_face = gauss_quad(0, 0, N + 1))
    VXY, EToV = uniform_mesh(rd.element_type, K)
    md = MeshData((VXY, EToV), rd; 
                is_periodic = true)

    equations = CompressibleEulerEquations2D(1.4)

    psi(u, normal, ::CompressibleEulerEquations2D) = u[2] * normal[1] + u[3] * normal[2]

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

    function mymul!(du::StructArray, A, u, alpha=1, beta=0)
        # return LinearAlgebra.mul!(du, A, u, alpha, beta)
        #mul_by! = (du, u) -> LinearAlgebra.mul!(du, A, u, alpha, beta)
        mul_by! = (du, u) -> Octavian.matmul!(du, A, u, alpha, beta)
        StructArrays.foreachfield(mul_by!, du, u)
        return du
    end

    regularized_ratio(a, b) = a * b / (b^2 + 1e-16)

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

    function compute_visc(u, params)

        (; rd, md, equations) = params
        (; invMQTr, invMQTs) = params
        (; K_visc) = params

        du = similar(u)

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
        uM = vf
        @batch for i in eachindex(uM) 
            uM[i] = entropy2cons(vf[i], equations)
        end
        uP = uM[md.mapP]
        if all(md.is_periodic) == false
            @. uP[md.mapB] = initial_condition(md.xf[md.mapB], md.yf[md.mapB], equations)
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
    
        # Calc piecewise constant visc dissipation
        unscaled_viscous_dissipation = zeros(md.num_elements)
        @batch for e in eachindex(unscaled_viscous_dissipation)
            val = zero(eltype(unscaled_viscous_dissipation))
            for i in axes(md.wJq, 1)
                val += md.wJq[i, e] * (dot(sigma_x[i, e], K_visc[i, e], sigma_x[i, e]) + 
                                       dot(sigma_y[i, e], K_visc[i, e], sigma_y[i, e]))
            end
            unscaled_viscous_dissipation[e] = val
        end

        # element-local K_visc                   
        elementwise_viscosity = vec(map((a,b) -> max(0.0, -min(0, a) * b / (1e-16 + b * b)), 
                                    entropy_errors, unscaled_viscous_dissipation))

        return entropy_errors, unscaled_viscous_dissipation, elementwise_viscosity
    end

    u = StructArray{SVector{nvariables(equations), Float64}}(ntuple(_ -> similar(md.x), 4))
    u .= rd.Pq * initial_condition_test.(md.xq, md.yq, 0, equations)    
    interface_flux = flux_lax_friedrichs

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
                sigmaM = similar(rd.Vf * u)
            )


    entropy_errors, unscaled_viscous_dissipation, elementwise_viscosity = 
        compute_visc(u, params)
    @show maximum(-min.(0, entropy_errors)), maximum(elementwise_viscosity)
    # @show unscaled_viscous_dissipation
    
    push!(max_elementwise_visc, maximum(elementwise_viscosity))
    push!(max_visc_dissipation, maximum(unscaled_viscous_dissipation))
    push!(min_visc_dissipation, minimum(unscaled_viscous_dissipation))
    push!(max_entropy_errors, maximum(abs.(entropy_errors)))
    # push!(max_entropy_errors, maximum(abs.(-min.(0, entropy_errors))))
end

using LaTeXStrings
d = 2
h = 0.5 .^ (1:length(max_entropy_errors))
plot(h, max_elementwise_visc, marker=:dot, linewidth=2, label=L"\epsilon_k(u_h)")
plot!(h, max_entropy_errors, marker=:square, linewidth=2,  label=L"\sigma_k(u_h)")

# smooth rates
r1 = (2 * N + 1)
C = 3 * max_elementwise_visc[end] / h[end]^r1
plot!(h, C * h .^ r1, linestyle=:dash, linewidth=2, label=L"O(h^{2N+1})")

r2 = (2 * N + 2 + d)
C = .3 * max_entropy_errors[end] / h[end]^r2
plot!(h, C * h .^ r2, linestyle=:dot, linewidth=2, label=L"O(h^{2N+4})")

# # discontinuous rates
# plot!(h, .05 * h, linestyle=:dash, linewidth=2, label=L"O(h)")

xlabel!("Mesh size")
plot!(xaxis=:log, yaxis=:log, 
      xguidefontsize=14, yguidefontsize=14, tickfontsize=14, 
      legend=:bottomright, legendfontsize=14, dpi=500)

