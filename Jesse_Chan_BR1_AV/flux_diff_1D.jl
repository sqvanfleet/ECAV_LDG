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
        # return prim2cons(SVector(0.125, 0.0, .1), equations)
        return prim2cons(SVector(0.0125, 0.0, .01), equations)
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
rd = RefElemData(Line(), N; quad_rule_vol=gauss_quad(0, 0, N+1))
# rd = RefElemData(Line(), SBP(), N)
(VX,), EToV = uniform_mesh(rd.element_type, 50)

# @. VX = 0.5 * (1 + VX) # modified Sod
# init_condition = initial_condition_modified_sod
# tspan = (0.0, 0.2) # modified Sod

@. VX = 5 * VX # for Shu Osher
init_condition = initial_condition_shu_osher
tspan = (0.0, 1.8) # Shu-Osher

md = MeshData(((VX,), EToV), rd)

# for density wave
# init_condition = initial_condition_density_wave
# tspan = (0.0, 10.0) # density wave
# md = MeshData(((VX,), EToV), rd; is_periodic=true) 

equations = CompressibleEulerEquations1D(1.4)

psi(u, ::CompressibleEulerEquations1D) = u[2]
psi(u, normal, ::CompressibleEulerEquations1D) = u[2] * normal

regularized_ratio(a, b; tol=1e-14) = a * b / (b^2 + tol)

using Trixi.ForwardDiff
dudv(v, equations) = ForwardDiff.jacobian(v -> entropy2cons(v, equations), v)

function rhs!(du, u, params, t)

    (; rd, md, equations, interface_flux) = params

    # compute entropy variables
    v = rd.Pq * cons2entropy.(rd.Vq * u, equations)

    uM = entropy2cons.(rd.Vf * v, equations)
    uP = uM[md.mapP]

    u_q = entropy2cons.(rd.Vq * v, equations)
    u_h = [u_q; uM]

    if md.is_periodic[1] == false
        (; init_condition) = params
        uP[1] = init_condition(md.x[1], equations)
        uP[end] = init_condition(md.x[end], equations)
    end

    # surface terms
    du .= rd.Vf' * (@. interface_flux(uM, uP, SVector(md.nx), equations))

    # volume terms
    Q = rd.Pq' * rd.M * rd.Dr * rd.Pq
    E = rd.Vf * rd.Pq 
    B = diagm([-1; 1])
    VhTr = [rd.Vq; rd.Vf]'
    Qh_skew = [Q - Q' E' * B; 
               -B*E zeros(2, 2)]
    du_vol = similar(u_h[:, 1])
    for e in axes(du, 2)
        fill!(du_vol, zero(eltype(du_vol)))
        for i in axes(u_h, 1), j in axes(u_h, 1)
            u_i = u_h[i, e]
            u_j = u_h[j, e]
            du_vol[i] += Qh_skew[i,j] * flux_ranocha(u_i, u_j, 1, equations)            
        end
        view(du, :, e) .+= VhTr * du_vol
    end
    du .= inv(rd.M) * (-du ./ md.J)
    return du    
end

u = rd.Pq * init_condition.(rd.Vq * md.x, equations)

interface_flux = flux_lax_friedrichs
# interface_flux = flux_ranocha
# interface_flux = flux_central

params = (; rd, md, equations, interface_flux, init_condition)
ode = ODEProblem(rhs!, u, tspan, params)
sol = solve(ode, SSPRK43(), 
            dt = 1e-8,
            abstol=1e-8, reltol=1e-6,
            saveat=LinRange(tspan..., 100), 
            callback=AliveCallback(alive_interval=100))

u = sol.u[end]
pad_nans(u) = vec([u; fill(NaN, 1, size(u, 2))])

# using LaTeXStrings
# plot(pad_nans(rd.Vp * md.x), pad_nans(rd.Vp * getindex.(u, 1)), linewidth=2,
#       label=L"N=%$N, K=%$(md.num_elements)")
# plot!(legend=:topright, dpi=500, xformatter=:none, yformatter=:none, legendfontsize=14)

using MAT
using LaTeXStrings
vars = matread("weno5_shuosher.mat")   
plot(vec(vars["x"]), vec(vars["rho"]), label="WENO", ylims=(0.5, 5.5), linewidth=2)
plot!(pad_nans(rd.Vp * md.x), pad_nans(rd.Vp * getindex.(u, 1)), linewidth=2,
      label=L"N=%$N, K=%$(md.num_elements)")
plot!(legend=:bottomleft, xformatter=:none, yformatter=:none, dpi=500)
plot!(dpi=500, xformatter=:none, yformatter=:none, legendfontsize=14)
# plot!(xlims=(-4, 3), ylims=(2.5, 5)) # zoom 

