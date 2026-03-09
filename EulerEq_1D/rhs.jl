
include("initial_conditions.jl")
include("computing_quanities.jl")

psi(u, normal, ::CompressibleEulerEquations1D) = u[2] * normal
dudv(v, equations) = ForwardDiff.jacobian(v -> entropy2cons(v, equations), v)

function ramp(s,s0,xi,epsilon_0)
    if s < s0 - xi
        return 0.0
    elseif abs(s-s0) < xi
        return epsilon_0/2 * (1 + sin(pi/(2*xi) * (s - s0)))
    else
        return epsilon_0
    end
end

function rhs_shock_cap!(du_voa, u_voa, params, t)
    du = parent(du_voa)
    u = parent(u_voa)

    (; rd, md, equations) = params

    push!(params.t, t)

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
    # vM = rd.Vf * v #nodal
    # vP = vM[md.mapP] #nodal
    vM = cons2entropy.(uM, equations) #modal
    vP = cons2entropy.(uP, equations) #modal
    interface_flux = @. params.interface_flux(uM, uP, SVector(md.nxJ), equations) * md.Jf
    duvol = md.rxJ .* (invMQTr * rd.Pq * flux.(rd.Vq * u, 1, equations)) 
    du .= duvol + rd.LIFT * interface_flux

    du_Euler_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * (du./md.J)).^2))
    push!(params.du_Euler_norm, du_Euler_norm)

    # use ρ as indicator
    rho = getindex.(u, 1)
    indicator = pressure.(u, equations) .* rho
    u_modal = rd.VDM \ indicator
    EN = vec(sum(u_modal.^2, dims=1))
    ENm1 = vec(sum(u_modal[1:end-1, :].^2, dims=1))
    S = max.(u_modal[end, :].^2 ./ EN, u_modal[end-1, :].^2 ./ ENm1)
    s = log10.(S);
    element_sizes = transpose(maximum(md.x, dims=1) .- minimum(md.x, dims=1))
    epsilon_0 = 0.5*element_sizes/rd.N
    s_0 = -(15/2)*log10(rd.N)
    xi = (7/2)*log10(rd.N)

    epsilon = @. ramp(s, s_0, xi, epsilon_0)

    if params.AV_type == :BR1
        beta = 0.0
    elseif params.AV_type == :LDG
        beta = @. sign(md.nxJ)
    end


    theta = (md.rxJ .* (invMQTr * v) + rd.LIFT * 
            (@. 0.5 * (vP + vM) * md.nxJ - 0.5 * beta * (vP - vM) * md.nxJ)) ./ md.J
    
    K = dudv.(cons2entropy.(rd.Vq * u, equations), equations)

    sigma = rd.Pq * (K .* (rd.Vq * theta))

    sigma = sigma * Diagonal(vec(epsilon))
    sigma_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * sigma).^2))
    
    sigmaM = rd.Vf * sigma
    sigmaP = sigmaM[md.mapP]

    du_visc = md.rxJ .* (invMQTr * sigma) + 
        rd.LIFT * (@. 0.5 * (sigmaP + sigmaM) * md.nxJ + 0.5 * beta * (sigmaP - sigmaM) * md.nxJ) ./ md.J
    
    du_visc_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * du_visc).^2))
    
    du .-= md.rxJ .* (invMQTr * sigma) + 
        rd.LIFT * (@. 0.5 * (sigmaP + sigmaM) * md.nxJ + 0.5 * beta * (sigmaP - sigmaM) * md.nxJ)
    @. du /= -md.J

    du_Euler_AV_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * du).^2))

    if params.initial_condition_type == :density_wave
        u_exact = my_initial_condition_density_wave.(md.xq .- 0.1*t, equations)
        error = rd.Vq * parent(u) - u_exact
        num = sqrt(sum(md.wJq .* norm.(error).^2))
        den = sqrt(sum(md.wJq .* norm.(u_exact).^2))
        L2_error = num
        push!(params.L2_error, L2_error)
        l2_epsilon = sqrt(sum(epsilon.^2))
        push!(params.l2_epsilon, l2_epsilon)
    end

    push!(params.du_Euler_AV_norm, du_Euler_AV_norm)
    push!(params.max_epsilon, maximum(epsilon))
    push!(params.du_visc_norm, du_visc_norm)
    push!(params.sigma_norm, sigma_norm)

    return du

end

function rhs!(du_voa, u_voa, params, t)
    
    du = parent(du_voa)
    u = parent(u_voa)

    (; rd, md, equations) = params

    push!(params.t, t)

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
    # vM = rd.Vf * v #nodal
    # vP = vM[md.mapP] #nodal
    vM = cons2entropy.(uM, equations) #modal
    vP = cons2entropy.(uP, equations) #modal
    interface_flux = @. params.interface_flux(uM, uP, SVector(md.nx), equations) * md.Jf
    duvol = md.rxJ .* (invMQTr * rd.Pq * flux.(rd.Vq * u, 1, equations)) 
    du .= duvol + rd.LIFT * interface_flux

    du_Euler_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * (du./md.J)).^2))
    push!(params.du_Euler_norm, du_Euler_norm)

    if params.AV_type == :BR1
        beta = 0.0
    elseif params.AV_type == :LDG
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
    #K = dudv.(rd.Vq * v_avg, equations)
    sigma = rd.Pq * (K .* (rd.Vq * theta))
    den = sum(md.wJq .* dot.(rd.Vq * sigma, rd.Vq * theta), dims=1)
    epsilon = @. num * den / (1e-14 + den^2)
    sigma = sigma * Diagonal(vec(epsilon))
    sigma_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * sigma).^2))
    
    sigmaM = rd.Vf * sigma
    sigmaP = sigmaM[md.mapP]

    du_visc = md.rxJ .* (invMQTr * sigma) + 
        rd.LIFT * (@. 0.5 * (sigmaP + sigmaM) * md.nxJ + 0.5 * beta * (sigmaP - sigmaM) * md.nxJ) ./ md.J
    
    du_visc_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * du_visc).^2))
    
    du .-= md.rxJ .* (invMQTr * sigma) + 
        rd.LIFT * (@. 0.5 * (sigmaP + sigmaM) * md.nxJ + 0.5 * beta * (sigmaP - sigmaM) * md.nxJ)
    @. du /= -md.J

    du_Euler_AV_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * du).^2))

    if params.initial_condition_type == :density_wave
        u_exact = my_initial_condition_density_wave.(md.xq .- 0.1*t, equations)
        error = rd.Vq * parent(u) - u_exact
        num = sqrt(sum(md.wJq .* norm.(error).^2))
        den = sqrt(sum(md.wJq .* norm.(u_exact).^2))
        L2_error = num
        push!(params.L2_error, L2_error)
        l2_epsilon = sqrt(sum(epsilon.^2))
        push!(params.l2_epsilon, l2_epsilon)
    end

    if params.initial_condition_type == :density_wave_convergence
        u_exact = density_wave_convergence.(md.xq .- 0.1*t, equations)
        error = rd.Vq * parent(u) - u_exact
        num = sqrt(sum(md.wJq .* norm.(error).^2))
        den = sqrt(sum(md.wJq .* norm.(u_exact).^2))
        L2_error = num/den
        push!(params.L2_error, L2_error)
    end

    push!(params.du_Euler_AV_norm, du_Euler_AV_norm)
    push!(params.max_epsilon, maximum(epsilon))
    push!(params.du_visc_norm, du_visc_norm)
    push!(params.sigma_norm, sigma_norm)

    return du
end

function rhs_DG_FEM!(du_voa, u_voa, params, t)
    
    du = parent(du_voa)
    u = parent(u_voa)

    (; rd, md, equations) = params

    push!(params.t, t)

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

    vM = cons2entropy.(uM, equations) 
    vP = cons2entropy.(uP, equations) 
    interface_flux = @. params.interface_flux(uM, uP, SVector(md.nx), equations) * md.Jf
    duvol = md.rxJ .* (invMQTr * rd.Pq * flux.(rd.Vq * u, 1, equations)) 
    du .= duvol + rd.LIFT * interface_flux

    du_Euler_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * (du./md.J)).^2))
    push!(params.du_Euler_norm, du_Euler_norm)

    @. du /= -md.J

    if params.initial_condition_type == :density_wave
        u_exact = my_initial_condition_density_wave.(md.xq .- 0.1*t, equations)
        error = rd.Vq * parent(u) - u_exact
        num = sqrt(sum(md.wJq .* norm.(error).^2))
        den = sqrt(sum(md.wJq .* norm.(u_exact).^2))
        L2_error = num
        push!(params.L2_error, L2_error)
    end

    if params.initial_condition_type == :density_wave_convergence
        u_exact = density_wave_convergence.(md.xq .- 0.1*t, equations)
        error = rd.Vq * parent(u) - u_exact
        num = sqrt(sum(md.wJq .* norm.(error).^2))
        den = sqrt(sum(md.wJq .* norm.(u_exact).^2))
        L2_error = num/den
        push!(params.L2_error, L2_error)
    end

    return du
end


function solve_ode(u, tspan, params, abstol = 1e-6, reltol = 1e-4, number_of_saves = 1000)
    # Solve the ODE problem
    ode = ODEProblem(rhs!, VectorOfArray(u), tspan, params)
    sol = solve(ode, SSPRK43(); abstol, reltol, 
                saveat=LinRange(tspan..., number_of_saves), 
                callback=AliveCallback(alive_interval=100))
    return sol
end