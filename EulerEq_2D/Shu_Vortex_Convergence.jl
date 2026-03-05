using StartUpDG
using LinearAlgebra 
using StaticArrays
using OrdinaryDiffEq, RecursiveArrayTools
using JLD2
using Trixi
using Trixi.ForwardDiff
using Plots

psi_1(u, normal, ::CompressibleEulerEquations2D) = u[2] * normal
psi_2(u, normal, ::CompressibleEulerEquations2D) = u[3] * normal

dudv(v, equations) = ForwardDiff.jacobian(v -> entropy2cons(v, equations), v)

function rhs!(du_voa, u_voa, params, t)
    
        du = parent(du_voa)
        u = parent(u_voa)

        (; rd, md, equations, AV_discretization) = params

        #compute entropy variables 
        v = rd.Pq * cons2entropy.(rd.Vq * u, equations)

        invMQTr = -rd.M \ (rd.Dr' * rd.M)
        invMQTs = -rd.M \ (rd.Ds' * rd.M)
        #uM = rd.Vf * u
        #The function entropy2cons is having trouble with some problems
        uM = entropy2cons.(rd.Vf * v, equations)
        uP = uM[md.mapP]
        vM = cons2entropy.(uM, equations)
        vP = cons2entropy.(uP, equations)

        interface_flux = @. params.interface_flux(uM, uP, SVector(md.nxJ, md.nyJ), equations)

        duvol = md.rxJ .* (invMQTr * rd.Pq * flux.(rd.Vq * u, 1, equations)) + md.sxJ .* (invMQTs * rd.Pq * flux.(rd.Vq * u, 1, equations)) +
                md.ryJ .* (invMQTr * rd.Pq * flux.(rd.Vq * u, 2, equations)) + md.syJ .* (invMQTs * rd.Pq * flux.(rd.Vq * u, 2, equations))

        du .=  duvol + rd.LIFT * interface_flux

        if AV_discretization == :LDG
                beta = @. sign(2*md.nxJ + md.nyJ)
        elseif AV_discretization == :BR1
                beta = 0.0
        end
        
        #compute thetas
        theta_1 = (md.rxJ.* (invMQTr * v) + md.sxJ.* (invMQTs * v) 
                        + rd.LIFT * (@. 0.5 * (vP + vM) * md.nxJ - 0.5 * beta * (vP - vM) * md.nxJ)) ./ md.J

        theta_2 = (md.ryJ.* (invMQTr * v) + md.syJ.* (invMQTs * v)
                        + rd.LIFT * (@. 0.5 * (vP + vM) * md.nyJ - 0.5 * beta * (vP - vM) * md.nyJ)) ./ md.J
        #Compute Epsilon
        term1 = dot.(v, rd.M * duvol)
        term2 = @. rd.wf * (psi_1(uM, md.nxJ, equations) + psi_2(uM, md.nyJ, equations))
        delta = sum(term1, dims=1) + sum(term2, dims=1)
        num = @. -min(0, delta)
        v_avg = repeat(0.5 * sum(Diagonal(rd.wq) * (rd.Vq * v), dims=1), rd.Np, 1)
        K = dudv.(rd.Vq * v_avg, equations)
        sigma_1 = rd.Pq * (K .* (rd.Vq * theta_1))
        sigma_2 = rd.Pq * (K .* (rd.Vq * theta_2))
        den = sum(md.wJq .* dot.(rd.Vq * sigma_1, rd.Vq * theta_1), dims=1) + 
              sum(md.wJq .* dot.(rd.Vq * sigma_2, rd.Vq * theta_2), dims=1)
        epsilon = @. num * den / (1e-14 + den^2)
        #Compute sigmas
        
        sigma_1 = sigma_1 * Diagonal(vec(epsilon))
        sigma_2 = sigma_2 * Diagonal(vec(epsilon))
        sigma_1M = rd.Vf * sigma_1
        sigma_1P = sigma_1M[md.mapP]
        sigma_2M = rd.Vf * sigma_2
        sigma_2P = sigma_2M[md.mapP]

    
        du .-= md.rxJ .* (invMQTr * sigma_1) + md.sxJ .* (invMQTs * sigma_1) +
               md.ryJ .* (invMQTr * sigma_2) + md.syJ .* (invMQTs * sigma_2) +
               rd.LIFT * (@. 0.5 * (sigma_1P + sigma_1M) * md.nxJ + 0.5 * beta * (sigma_1P - sigma_1M) * md.nyJ) +
               rd.LIFT * (@. 0.5 * (sigma_2P + sigma_2M) * md.nyJ + 0.5 * beta * (sigma_2P - sigma_2M) * md.nyJ)  

        @. du /= -md.J
    
        return u

end

function shu_isentropic_euler_vortex(x,y,t,L,gamma,equations::CompressibleEulerEquations2D)
    alpha = pi/4
    M_inf = sqrt(2.0/gamma)
    sigma = 1.0
    R = 1.0
    Beta = M_inf*5.0*sqrt(2.0)*exp(0.5)/(4*pi)
    ux = M_inf * cos(alpha)
    vx = M_inf * sin(alpha)
    #t = periodic_wrap(t, 0.0, 2.0*L*sqrt(gamma))
    x = periodic_wrap(x, -L, 2.0*L)
    y = periodic_wrap(y, -L, 2.0*L)
    r = sqrt(((x/R)-ux*t)^2 + ((y/R)-vx*t)^2)
    Omega(r) = Beta * exp(-r^2/(2.0*sigma^2)) 
    rho = (1.0 - 0.5 * (gamma - 1.0) * Omega(r)^2)^(1.0/(gamma - 1.0))
    u = M_inf * cos(alpha) - (y/R - vx*t)*Omega(r)
    v = M_inf * sin(alpha) + (x/R - ux*t)*Omega(r)
    p = (1/gamma)*(1.0 - 0.5 * (gamma - 1.0) * Omega(r)^2)^(gamma/(gamma - 1.0))
    return prim2cons(SVector(rho, u, v, p), equations)
end

gamma = 1.4
equations = CompressibleEulerEquations2D(gamma)

#Choose Modal or Nodal
#DG_type = :nodal
DG_type = :modal

#BR1 versus LDG switch
#AV_discretization = :BR1
AV_discretization = :LDG

function periodic_wrap(x, a, T)
    return mod(x - a, T) + a
end

function compute_L2_error(N, K1D)
    if DG_type == :nodal
        rd = RefElemData(Tri(), SBP(), N)
    elseif DG_type == :modal
        rd = RefElemData(Tri(), N)
    end

    VXY, EToV = uniform_mesh(rd.element_type, K1D)
    md = MeshData((VXY,EToV), rd; is_periodic = true)

    VXY, EToV = uniform_mesh(rd.element_type, K1D)
    L = 10.0 #Length of the domain
    VXY = VXY.*L #Scale the mesh
    md = MeshData((VXY,EToV), rd; is_periodic = true)
    x,y = md.x, md.y
    u = rd.Pq * shu_isentropic_vortex.(md.xq,md.yq,0.0,L,gamma,equations)

#     Final_time = 2.0*sqrt(gamma)*L - 5.0
    Final_time = 15.0
    tspan = (0.0, Final_time)
    du = similar(VectorOfArray(u))
    params = (; rd, md, equations, interface_flux = flux_lax_friedrichs, AV_discretization)
    ode = ODEProblem(rhs!, VectorOfArray(u), tspan, params)
    sol = solve(ode, SSPRK43(); abstol = 1e-8, reltol = 1e-6, 
    saveat=LinRange(tspan..., 50), 
    callback=AliveCallback(alive_interval=100))

    u = parent(sol.u[end])
    alpha = pi/4
    M_inf = sqrt(2.0/gamma)
    #Compute the exact solution
    ux = M_inf * cos(alpha)
    uy = M_inf * sin(alpha)
    u_exact = @. shu_isentropic_euler_vortex(md.xq-ux*Final_time,md.yq-uy*Final_time,0.0,L,gamma,equations)
    L2_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * u - u_exact).^2))
    L2_rel_norm = L2_norm / sqrt(sum(md.wJq .* norm.(rd.Vq * u_exact).^2))
    @show L2_norm
    @show L2_rel_norm   

    if DG_type == :nodal && AV_discretization == :LDG
            @save "Vortex_Data_10/N$(N)_K$(K1D)_T$(Final_time)_nodal_LDG.jld2" md rd sol u_exact L2_norm L2_rel_norm        
    elseif DG_type == :nodal && AV_discretization == :BR1
            @save "Vortex_Data_10/N$(N)_K$(K1D)_T$(Final_time)_nodal_BR1.jld2" md rd sol u_exact L2_norm L2_rel_norm
    elseif DG_type == :modal && AV_discretization == :BR1
            @save "Vortex_Data_10/N$(N)_K$(K1D)_T$(Final_time)_modal_BR1.jld2" md rd sol u_exact L2_norm L2_rel_norm
    elseif DG_type == :modal && AV_discretization == :LDG           
            @save "Vortex_Data_10/N$(N)_K$(K1D)_T$(Final_time)_modal_LDG.jld2" md rd sol u_exact L2_norm L2_rel_norm
    end

end

L2_errors = [compute_L2_error(N, K1D) for N in [1], K1D in [16, 32, 64]]