using StartUpDG
using LinearAlgebra: dot
using OrdinaryDiffEq, RecursiveArrayTools
using JLD2
using Trixi
using LaTeXStrings

#equations = LinearScalarAdvectionEquation1D(1.0)
equations = InviscidBurgersEquation1D()

N = 4
K1D = 32 
rd = RefElemData(Line(),SBP(), N)
md = MeshData(uniform_mesh(Line(), K1D), rd;
              is_periodic=true)
              
psi(u, normal, ::LinearScalarAdvectionEquation1D) = u[1]^2 / 2 * normal
psi(u, normal, ::InviscidBurgersEquation1D) = u[1]^3 / 6 * normal
function Trixi.flux_ec(u_ll, u_rr, normal::AbstractVector, 
                       ::InviscidBurgersEquation1D) 
    return flux_ec(u_ll, u_rr, 1, equations) * normal[1] 
    #return flux_lax_friedrichs(u_ll, u_rr, normal[1], equations)
end

function rhs!(du_voa, u_voa, params, t)

    du = parent(du_voa)
    u = parent(u_voa)

    (; rd, md, equations, AV_type) = params
    invMQTr = -rd.M \ (rd.Dr' * rd.M) 
    uM = rd.Vf * u 
    uP = uM[md.mapP] 

    interface_flux = @. params.interface_flux(uM, uP, SVector(md.nxJ), equations)
    duvol = md.rxJ .* (invMQTr * flux.(u,1,equations))
    du .= duvol + rd.LIFT * interface_flux # Right hand side for the, LIFT matrix has Minv built into it.

    if AV_type == :LDG
            beta = @. sign(md.nxJ)
    elseif AV_type == :BR1
            beta = 0.0
    end

    theta = (md.rxJ .* (invMQTr * u) + rd.LIFT * (@. 0.5 * (uP + uM) * md.nxJ - 0.5 * beta * (uP - uM) * md.nxJ)) ./ md.J
    
    #Compute Epsilon
    term1 = dot.(u, rd.M * duvol)
    term2 = @. rd.wf * (psi(uM,md.nxJ,equations))
    delta = sum(term1, dims = 1) + sum(term2, dims=2)
    num = @. -min(0, delta)
    sigma = rd.Pq * (rd.Vq * theta)
    den = sum(md.wJq .* dot.(rd.Vq * sigma, rd.Vq * theta), dims = 1)
    epsilon = @. num*den / (1e-14 + den^2)

    sigma = sigma * Diagonal(vec(epsilon))
    sigmaM = rd.Vf * sigma
    sigmaP = sigmaM[md.mapP]
    du .-= md.rxJ .* (invMQTr * sigma) + rd.LIFT * (@. 0.5 * (sigmaP + sigmaM) * md.nxJ + 0.5 * beta * (sigmaP - sigmaM) * md.nxJ)
    @. du /= -md.J
    
    return du
end

AV_type = :LDG
x = md.x
u = @. SVector(0.5 - sin(pi * x))
params = (; rd, md, equations, interface_flux = flux_lax_friedrichs, AV_type) #could also use interface_flux = flux_ec
tspan = (0, 1.0)
ode = ODEProblem(rhs!, VectorOfArray(u), tspan, params)
sol = solve(ode, SSPRK43(); abstol = 1e-6, reltol = 1e-4, 
            saveat=LinRange(tspan..., 50))

function compute_quantities(u,params)
    (; rd, md, equations, AV_type) = params
    invMQTr = -rd.M \ (rd.Dr' * rd.M) 
    uM = rd.Vf * u 
    uP = uM[md.mapP] 

    interface_flux = @. params.interface_flux(uM, uP, SVector(md.nxJ), equations)
    du = similar(u)
    du .= md.rxJ .* (invMQTr * flux.(u, 1, equations)) + rd.LIFT * interface_flux # Right hand side for the, LIFT matrix has Minv built into it.

    if AV_type == :LDG
            beta = @. sign(md.nxJ)
    elseif AV_type == :BR1
            beta = 0.0
    end

    #sigma = (md.rxJ .* (invMQTr * u) + rd.LIFT * (@. 0.5 * (uP + uM))) ./ md.J
    sigma = (md.rxJ .* (invMQTr * u) + rd.LIFT * (@. 0.5 * (uP + uM) * md.nxJ - 0.5 * beta * (uP - uM) * md.nxJ)) ./ md.J
    
    K1D = md.num_elements
    delta = zeros(K1D)
    num = zeros(K1D)
    den = zeros(K1D)
    epsilon = zeros(K1D)
    dvdxJ = md.rxJ .* (rd.Dr * u)
    term1 = @. -rd.wq * dot(dvdxJ, flux(u, 1, equations))
    term2 = @. rd.wf * psi(uM, md.nxJ, equations)
    for i in 1:K1D
        delta[i] = sum(term1[:,i]) + sum(term2[:,i])
        # Comput the numerator of eq (22)
        num[i] = -min(0, delta[i])
        # Compute the denomenator of eq (22)
        #den[i] = sum(rd.wq .* sigma[:,i].^2 .* md.J[:,i])
        den[i] = sum(rd.wq .* dot.(sigma[:,i], sigma[:,i]) .* md.J[:,i])
        epsilon[i] = num[i] * den[i] / (1e-14 + den[i]^2)
        sigma[:,i] .*= epsilon[i] 
    end

    sigmaM = rd.Vf * sigma
    sigmaP = sigmaM[md.mapP]
    du .-= md.rxJ .* (invMQTr * sigma) + rd.LIFT * (@. 0.5 * (sigmaP + sigmaM) * md.nxJ + 0.5 * beta * (sigmaP - sigmaM) * md.nxJ)
    @. du /= -md.J


    dSdt = sum(dot.(u,rd.M*(md.J .* du)))
    #dSdt = sum(u.*(rd.M*(md.J .* du)))
    return dSdt, epsilon, delta
end

using Plots
u = parent(sol.u[end])

epsilon_evolution_LDG = zeros(length(sol.t))
dSdt_evolution_LDG  = zeros(length(sol.t))
delta_evolution_LDG = zeros(length(sol.t))

for (i, u) in enumerate(sol.u)
    temp =  parent(u)# Ensure u is evaluated
    display(i)
    #Compute Epsilon
    dSdt, epsilon, delta = compute_quantities(temp, params)
    epsilon_evolution_LDG[i] = maximum(epsilon)
    delta_evolution_LDG[i] = maximum(delta)
    dSdt_evolution_LDG[i] = dSdt
end

p1 = plot(sol.t, dSdt_evolution_LDG,legend = false, xlabel = L"$t$", ylabel = L"$\frac{\mathrm{d}S}{\mathrm{d}{t}}$"
    ,yguidefontrotation = -90, left_margin = 10Plots.mm, title = L"Evolution of $\frac{\mathrm{d}S}{\mathrm{d}{t}}$")
display(p1)

p1 = plot(sol.t, epsilon_evolution_LDG, legend = false, xlabel = L"$t$", ylabel = L"$\varepsilon$"
    ,yguidefontrotation = -90, left_margin = 10Plots.mm, title = L"Evolution of $\varepsilon$")
display(p1)

p1 = plot(rd.Vp * md.x, rd.Vp * getindex.(u,1),yguidefontrotation = -90, left_margin = 10Plots.mm, legend = false, title = "Solution at the Final time",
        xlabel = "x", ylabel = "y")
display(p1)


