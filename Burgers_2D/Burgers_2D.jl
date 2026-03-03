using StartUpDG
using OrdinaryDiffEqSSPRK
using Plots
using Trixi: AliveCallback
using JLD2


f1(u) = u^2 / 2
f2(u) = u^2 / 2
psi1(u) = u^3 / 6
psi2(u) = u^3 / 6

function rhs!(du, u, params, t)
    (; rd, md, AV_type ) = params
    invMQTr = -rd.M \ (rd.Dr' * rd.M * rd.Pq)
    invMQTs = -rd.M \ (rd.Ds' * rd.M * rd.Pq)

    push!(params.t, t)

    uq = rd.Vq * u
    uM = rd.Vf * u
    uP = uM[md.mapP]
    if flux_type == :EC
        interface_flux = @. 1/6 * (uP^2 + uM * uP + uM^2) * (md.nxJ + md.nyJ)
    end

    if flux_type == :LxW
        lambda = @. max(abs(uP), abs(uM))
        interface_flux = @. 0.5 * (f1(uP) + f1(uM)) * md.nxJ + 
                            0.5 * (f2(uP) + f2(uM)) * md.nyJ +
                        - 0.5 * lambda * (uP - uM) * md.Jf
    end 


    du_vol = md.rxJ .* (invMQTr * f1.(uq)) + md.sxJ .* (invMQTs * f1.(uq)) +
             md.ryJ .* (invMQTr * f2.(uq)) + md.syJ .* (invMQTs * f2.(uq))
    du .= du_vol + rd.LIFT * interface_flux
    
    if AV_type == :LDG
        beta = @. sign(2*md.nx + md.ny)
    elseif AV_type == :BR1
        beta = 0
    end

    sigma_x = (md.rxJ .* (invMQTr * uq) + md.sxJ .* (invMQTs * uq) + 
               rd.LIFT * (@. (0.5 * (uP + uM) - 0.5 * beta * (uP - uM)) * md.nxJ)) ./ md.J
    sigma_y = (md.ryJ .* (invMQTr * uq) + md.syJ .* (invMQTs * uq) + 
               rd.LIFT * (@. (0.5 * (uP + uM) - 0.5 * beta * (uP - uM)) * md.nyJ)) ./ md.J

    # calc AV coefficient
    delta = sum(u .* (rd.M * du_vol), dims=1) + 
            sum((@. rd.wf * (psi1(uM) * md.nxJ + psi2(uM) * md.nyJ)), dims = 1)
    num = @. -min(0, delta)
    den = sum(md.wJq .* ((rd.Vq * sigma_x) .^ 2 + (rd.Vq * sigma_y) .^ 2), dims = 1)
    epsilon = @. num * den / (eps(den) + den^2)
    #epsilon = @. num * den / (1e-14 + den^2)

    
    sigma_x .*= epsilon
    sigma_y .*= epsilon

    invMQTr = invMQTr * rd.Vq
    invMQTs = invMQTs * rd.Vq
    sigma_xM = rd.Vf * sigma_x
    sigma_yM = rd.Vf * sigma_y
    sigma_xP = sigma_xM[md.mapP]
    sigma_yP = sigma_yM[md.mapP]
    interface_flux = @. (0.5 * (sigma_xP + sigma_xM) + 0.5 * beta * (sigma_xP - sigma_xM)) * md.nxJ + 
                        (0.5 * (sigma_yP + sigma_yM) + 0.5 * beta * (sigma_yP - sigma_yM)) * md.nyJ
    du .-= (md.rxJ .* (invMQTr * sigma_x) + md.sxJ .* (invMQTs * sigma_x)
          + md.ryJ .* (invMQTr * sigma_y) + md.syJ .* (invMQTs * sigma_y)
          + rd.LIFT * interface_flux)

    push!(params.max_epsilon, maximum(epsilon))
    
    @. du /= -md.J

    # @show extrema(epsilon)

    dSdt = sum(u .* (rd.M * (md.J .* du)))
    if dSdt > 100 * eps()
        @show dSdt
    end
    push!(params.dSdt, dSdt)
end

N = 3
K1D = 30
flux_type = :EC
DG_type = :modal  # Set to modal as requested

# Define the cases we want to run
av_cases = [:LDG, :BR1]

for current_AV in av_cases
    let rd, md, u, params, tspan, ode, sol, save_path 
        println("\n" * "="^40)
        println("Running Case: DG_type = $DG_type | AV_type = $current_AV")
        println("="^40)

        rd = RefElemData(Tri(), N) 
        md = MeshData(uniform_mesh(Tri(), K1D), rd; is_periodic = true)

        u = rd.Pq * exp.(-25*(md.xq.^2 + md.yq.^2)) 

        params = (; rd, md, AV_type = current_AV, DG_type, 
                   t = Float64[], max_epsilon = Float64[], dSdt = Float64[])
        
        tspan = (0.0, 0.77)
        ode = ODEProblem(rhs!, u, tspan, params)

        sol = solve(ode, SSPRK43(); abstol = 1e-6, reltol = 1e-4, 
                    saveat = LinRange(tspan..., 500),
                    callback = AliveCallback(alive_interval = 100))

        file_name = "N$(N)_K$(K1D)_$(current_AV)_$(DG_type)_$(flux_type).jld2"
        save_path = joinpath(@__DIR__, "Data", file_name)
        
        mkpath(dirname(save_path))
        
        println("Saving results to: $save_path")
        @save save_path md rd sol
    end
end

