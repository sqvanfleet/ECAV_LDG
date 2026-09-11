using OrdinaryDiffEq
using StartUpDG
using Plots, LinearAlgebra
using Trixi
using StaticArrays
using Polyester, Octavian
using StructArrays

N = 3
cells_per_dimension = 64
rd = RefElemData(Tri(), N)
VXY, EToV = uniform_mesh(rd.element_type, cells_per_dimension)
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

regularized_ratio(a, b) = a * b / (b^2 + 1e-14)

function rhs!(du, u, params, t)

    (; rd, md, equations, interface_flux) = params
    (; invMQTr, invMQTs) = params
    
    # interpolate to quad points, project entropy variables
    uq = rd.Vq * u
    vq = similar(uq)
    @batch for i in eachindex(vq)
        vq[i] = cons2entropy(uq[i], equations)
    end
    v = rd.Pq * vq

    # interpolate
    vf = rd.Vf * v   
    uM = vf
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

    (; correction_version) = params    
    if correction_version == :Abgrall
        u_avg = params.eTrM * u
        v_avg = params.eTrM * v
        v_osc = rd.Vq * v .- v_avg
        v_osc_local = similar(v_osc[:,1])
        @batch for e in axes(u, 2)
            K = dudv_cons_vars(u_avg[e], equations)
            # K = dudv_cons_vars(entropy2cons(v_avg[e], equations), equations)
            a = -min(0, entropy_errors[e])
            v_osc_local .= view(v_osc, :, e)
            b = 0.0
            for i in eachindex(rd.wq)
                b += rd.wq[i] * v_osc_local[i]' * K * v_osc_local[i]
            end
            mymul!(view(du, :, e), rd.Pq, 
                   regularized_ratio(a, b) * v_osc_local, 1, 1)
        end
    else # if correction == :Offner
        u_avg = params.eTrM * u
        VqDr, VqDs = rd.Vq * rd.Dr, rd.Vq * rd.Ds
        dvdr = VqDr * v
        dvds = VqDs * v
        GtGvr, GtGvs = ntuple(_ -> similar(dvdr[:, 1]), 2)
        @batch for e in axes(u, 2)
            K = dudv_cons_vars(u_avg[e], equations)
            a = -min(0, entropy_errors[e])
            
            G = SMatrix{2, 2}(md.rxJ[1, e], md.sxJ[1, e], 
                              md.ryJ[1, e], md.syJ[1, e]) / md.J[1, e]
            GtG = G' * G 

            for i in eachindex(GtGvr)
                GtGvr[i] = GtG[1, 1] * dvdr[i, e] + GtG[1, 2] * dvds[i, e]
                GtGvs[i] = GtG[2, 1] * dvdr[i, e] + GtG[2, 2] * dvds[i, e]
            end
            b = 0.0
            for i in eachindex(rd.wq)
                b += md.J[1, e] * rd.wq[i] * (GtGvr[i]' * K * GtGvr[i] + GtGvs[i]' * K * GtGvs[i])
            end
            mymul!(view(du, :, e), params.invM * (VqDr'), 
                   regularized_ratio(a, b) * md.J[1, e] * rd.wq .* GtGvr, 1, 1)
            mymul!(view(du, :, e), params.invM * (VqDs'), 
                   regularized_ratio(a, b) * md.J[1, e] * rd.wq .* GtGvs, 1, 1)                   
        end
    end

    # now add transport surface terms
    (; convective_interface_flux) = params
    @batch for i in eachindex(convective_interface_flux)
        convective_interface_flux[i] = 
            interface_flux(uM[i], uP[i], SVector(md.nx[i], md.ny[i]), equations) * md.Jf[i]
    end
    mymul!(du, rd.LIFT, convective_interface_flux, 1, 1)  

    @batch for i in eachindex(du)
        du[i] /= -md.J[i]
    end
    
end

# Riemann problem
function initial_condition_Riemann_problem(x, y, t, equations::CompressibleEulerEquations2D)

    if (0 < x < 1) & (0 < y < 1)
        rho = .5313
        u, v = 0, 0
        p = 0.4
    elseif (-1 < x < 0) & (0 < y < 1)
        rho = 1
        u, v = .7276, 0
        p = 1
    elseif (-1 < x < 0) & (-1 < y < 0)
        rho = .8
        u, v = 0, 0
        p = 1
    elseif (0 < x < 1) & (-1 < y < 0)        
        rho = 1
        u, v = 0, .7276
        p = 1
    else 
        @show x, y
    end
    return prim2cons(SVector(rho, u, v, p), equations)
end

init_cond = initial_condition_Riemann_problem

u = StructArray{SVector{nvariables(equations), Float64}}(ntuple(_ -> similar(md.x), 4))
u .= rd.Pq * init_cond.(md.xq, md.yq, 0.0, equations)

interface_flux = flux_lax_friedrichs

params = (; rd, md, equations, interface_flux, 
            v=similar(u), 
            gradients = (similar(u), similar(u)),
            viscosity=ones(size(u, 2)), 
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
            eTrM = 0.5 * sum(rd.M, dims=1), 
            integrand = zeros(size(md.xq, 1)),
            invM = inv(rd.M),
            correction_version = :Abgrall,
            # correction_version = :Offner,
         )

tspan = (0.0, .25)

ode = ODEProblem(rhs!, u, tspan, params)
sol = solve(ode, SSPRK43(), 
            dt = 1e-8,
            abstol=1e-6, reltol=1e-4,
            saveat=LinRange(tspan..., 50), 
            callback=AliveCallback(alive_interval=100))

u = sol.u[end]
scatter(vec(rd.Vp * md.x), vec(rd.Vp * md.y), 
        zcolor=vec(rd.Vp * getindex.(u, 1)), 
        legend=false, ms=.5, msw=0, ratio=1)

interp = vandermonde(rd.element_type, rd.N, equi_nodes(rd.element_type, rd.N)...) / rd.VDM
pdata = [interp * getindex.(u, 1), interp * pressure.(u, equations)]
# pdata = [interp * getindex.(u, 1), interp * pressure.(u, equations), interp * repeat(params.viscosity', rd.Np, 1)]
num_elems = Int(sqrt(md.num_elements ÷ 2))
# vtu_name = MeshData_to_vtk(md, rd, pdata, ["rho", "p"], "Riemann_N$(N)_K$(num_elems)_Offner", true)

