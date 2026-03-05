using StartUpDG
using LinearAlgebra 
using StaticArrays
using OrdinaryDiffEq, RecursiveArrayTools
using JLD2
using Trixi
using Trixi.ForwardDiff
using Plots

include("initial_conditions.jl")

psi_1(u, normal, ::CompressibleEulerEquations2D) = u[2] * normal
psi_2(u, normal, ::CompressibleEulerEquations2D) = u[3] * normal

dudv(v, equations) = ForwardDiff.jacobian(v -> entropy2cons(v, equations), v)

function rhs!(du_voa, u_voa, params, t)
    
        du = parent(du_voa)
        u = parent(u_voa)

        (; rd, md, equations, AV_discretization, gamma) = params

        #Get u at the quadrature points
        uq = rd.Vq * u

        #compute entropy variables 
        v = rd.Pq * cons2entropy.(uq, equations)

        invMQTr = -rd.M \ (rd.Dr' * rd.M)
        invMQTs = -rd.M \ (rd.Ds' * rd.M)
        uM = entropy2cons.(rd.Vf * v, equations)
        uP = uM[md.mapP]

        if all(md.is_periodic .== (true,false)) 

                tol = 100 * eps()
                top_bottom_boundary = findall(@. abs(md.yf[md.mapB] - 1.0) < tol || abs(md.yf[md.mapB]) < tol)
                # right_boundary = findall(@. abs(md.xf[md.mapB] - 2.0) < tol)
                # left_boundary = findall(@. abs(md.xf[md.mapB]) < tol)
                # left_right_boundary = findall(@. abs(md.xf[md.mapB] - 2.0) < tol || abs(md.xf[md.mapB]) < tol)


                # impose wall BC
                wall_indices = top_bottom_boundary
                for i in eachindex(wall_indices)
                        rho, v1, v2, p = cons2prim(uM[md.mapB[i]], equations)
                        nx, ny = md.nx[md.mapB[i]], md.ny[md.mapB[i]] 
                        v_normal = v1 * nx + v2 * ny
                        uP[md.mapB[i]] = prim2cons(SVector(rho, v1 , v2 - 2*v_normal*ny, p), equations)
                end

                # wall_indices_lr = left_right_boundary
                # for i in eachindex(wall_indices_lr)
                #         rho, v1, v2, p = cons2prim(uM[md.mapB[i]], equations)
                #         nx, ny = md.nx[md.mapB[i]], md.ny[md.mapB[i]] 
                #         v_normal = v1 * nx + v2 * ny
                #         uP[md.mapB[i]] = prim2cons(SVector(rho, v1 - 2*v_normal*nx, v2, p), equations)
                # end

                # # impose free-stream BCs on the left boundary
                # free_stream = left_boundary
                # for i in eachindex(free_stream)
                #         uP[md.mapB[i]] = Shock_vortex_interaction(md.xf[md.mapB[i]], md.yf[md.mapB[i]], equations, gamma)
                # end

                # # impose subsonic outflow BCs on the right boundary
                # subsonic_outflow = right_boundary
                # for i in eachindex(subsonic_outflow)
                #         rho_init, v1_init, v2_init, p_init = 
                #                 cons2prim(Shock_vortex_interaction(md.xf[md.mapB[i]], md.yf[md.mapB[i]], equations, gamma),equations) 
                #         rho, v1, v2, p = cons2prim(uM[md.mapB[i]], equations)
                #         uP[md.mapB[i]] = prim2cons(SVector(rho, v1, v2, p_init), equations)
                # end                
                        
        end

        vM = cons2entropy.(uM, equations)
        vP = cons2entropy.(uP, equations)

        interface_flux = @. params.interface_flux(uM, uP, SVector(md.nx, md.ny), equations)*md.Jf


        duvol = md.rxJ .* (invMQTr * rd.Pq * flux.(uq, 1, equations)) +
                 md.sxJ .* (invMQTs * rd.Pq * flux.(uq, 1, equations)) +
                 md.ryJ .* (invMQTr * rd.Pq * flux.(uq, 2, equations)) +
                 md.syJ .* (invMQTs * rd.Pq * flux.(uq, 2, equations))

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
        K = dudv.(cons2entropy.(uq, equations), equations)
        # K = dudv.(rd.Vq * v_avg, equations)
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
               rd.LIFT * (@. 0.5 * (sigma_1P + sigma_1M) * md.nxJ + 0.5 * beta * (sigma_1P - sigma_1M) * md.nxJ) +
               rd.LIFT * (@. 0.5 * (sigma_2P + sigma_2M) * md.nyJ + 0.5 * beta * (sigma_2P - sigma_2M) * md.nyJ)  

        @. du /= -md.J
    
        return du

end

gamma = 1.4
equations = CompressibleEulerEquations2D(gamma)

#Choose the 2D Euler equations problem
# initial_condition_type = :density_wave
# initial_condition_type = :shu_isentropic_vortex
# initial_condition_type = :Riemann_problem_config_1
initial_condition_type = :Shock_vortex_interaction

#Choose Modal or Nodal
#DG_type = :nodal
DG_type = :modal

#Mesh Parameters
N = 2 #Order of Polynomial
K1D = 16 #Number of elements in each direction

#BR1 versus LDG switch
#AV_discretization = :BR1
AV_discretization = :LDG

if DG_type == :nodal
        rd = RefElemData(Tri(), SBP(), N)
elseif DG_type == :modal
        rd = RefElemData(Tri(), N)
end

if initial_condition_type == :Shock_vortex_interaction
    (VX, VY), EToV = uniform_mesh(rd.element_type, 2*K1D, K1D)
    VX = VX .+ 1.0
    VY = 0.5 * (VY .+ 1.0)
    md = MeshData((VX, VY), EToV, rd; is_periodic = (true, false))

    u_init = rd.Pq * Shock_vortex_interaction.(md.xq, md.yq, equations, gamma)

    Final_time = 0.7
    tspan = (0.0, Final_time) 
    params = (; rd, md, equations, interface_flux = flux_hllc, AV_discretization, gamma)
    
    ode = ODEProblem(rhs!, VectorOfArray(u_init), tspan, params)
    sol = solve(ode, SSPRK43(); abstol = 1e-6, reltol = 1e-4, 
                saveat = LinRange(tspan..., 50), 
                callback = AliveCallback(alive_interval=100))

    data_path = joinpath(@__DIR__, "Shock_vortex_interaction_Data")
    mkpath(data_path) 
    xp, yp = md.x, md.y

    # Define the indices we want to plot: Start, Middle-ish, and End
    plot_indices = [1, length(sol.t)-35, length(sol.t)]

    for i in plot_indices
        u_curr = parent(sol.u[i])
        t_val = sol.t[i]
        rho = getindex.(u_curr, 1)

        # --- 1. Calculate & Plot Schlieren ---
        drhodx = (md.rxJ .* (rd.Dr * rho) + md.sxJ .* (rd.Ds * rho)) ./ md.J
        drhody = (md.ryJ .* (rd.Dr * rho) + md.syJ .* (rd.Ds * rho)) ./ md.J
        g = sqrt.(drhodx.^2 + drhody.^2)
        
        gmax, gmin = maximum(g), minimum(g)
        rho_schl = @. exp(-10 * (g - gmin) / (gmax - gmin))
        plots_path = joinpath(@__DIR__, "Shock_vortex_interaction_plots")
        mkpath(plots_path)
        p_schl = scatter(vec(xp), vec(yp), zcolor=vec(rho_schl),
                         msw=0, ms=0.5, legend=false, ratio=1, cam=(0,90), 
                         c=:red, colorbar=true, guidefontsize=20, tickfontsize=20)
        savefig(p_schl, joinpath(plots_path, "schlieren_t_$(t_val).png"))
        export_to_vtk(rd, md, [rho_schl], joinpath(plots_path, "schlieren_t_$(t_val).vtu"))
        

        # --- 2. Plot Density ---
        p_dens = scatter(vec(xp), vec(yp), zcolor=vec(rho), clims=(0.85, 1.20), 
                         msw=0, ms=0.5, legend=false, ratio=1, cam=(0,90),
                         c=:red, colorbar=true, guidefontsize=20, tickfontsize=20)
        savefig(p_dens, joinpath(plots_path, "density_t_$(t_val).png"))
        export_to_vtk(rd, md, [rho], joinpath(plots_path, "density_t_$(t_val).vtu"))
    end

    # --- 3. Save Data ---
    filename = joinpath(data_path, "N$(N)_K$(2*K1D)x$(K1D)_T$(Final_time)_$(DG_type)_$(AV_discretization).jld2")
    @save filename md rd sol
end


#Density wave problem
if initial_condition_type == :density_wave
        N_list = [1]
        K_list = [4,8,16,32,64]


        DG_type = :modal
        AV_discretization = :BR1

        for N in N_list
                for K1D in K_list
                        #Mesh Parameters
                        # N = 1 #Order of Polynomial
                        # K1D = 8 #Number of elements in each direction
                        display("N = $(N), K = $(K1D)")
                        if DG_type == :nodal
                                local rd = RefElemData(Tri(), SBP(), N)
                        elseif DG_type == :modal
                                local rd = RefElemData(Tri(), N)
                        end

                        local VXY, EToV = uniform_mesh(rd.element_type, K1D)
                        local md = MeshData((VXY,EToV), rd; is_periodic = true)
                        local u = rd.Pq * density_wave.(md.xq,md.yq,0.0,equations)

                        # display(scatter(vec(rd.Vp * md.x), vec(rd.Vp * md.y), zcolor=vec(rd.Vp * getindex.(u,1)), msw=0, ms = 1, legend=false, ratio=1, cam=(0,90),
                        #         title="DG_Density", xlabel="x", ylabel="y", c=:red))

                        local Final_time = 1.7
                        local tspan = (0.0, Final_time)
                        local du = similar(VectorOfArray(u))
                        local params = (; rd, md, equations, interface_flux = flux_lax_friedrichs, AV_discretization, gamma)
                        local ode = ODEProblem(rhs!, VectorOfArray(u), tspan, params)
                        local sol = solve(ode, SSPRK43(); abstol = 1e-8, reltol = 1e-6, 
                                        saveat=LinRange(tspan..., 50), 
                                        callback=AliveCallback(alive_interval=100))
                        local u = parent(sol.u[end])
                        local xp, yp, up = rd.Vp * md.x, rd.Vp * md.y, rd.Vp * getindex.(u,1) 
                        
                        # # higher accuracy quadrature rule
                        local rq, sq, wq = quad_nodes(rd.element_type, rd.N)
                        local Vq = vandermonde(rd.element_type, rd.N, rq, sq) / rd.VDM
                        local xq, yq = Vq * md.x, Vq * md.y
                        # local u_exact = density_wave.(xq, yq, sol.t[end], equations)
                        # local error = Vq * parent(u) - u_exact

                        # Defult quadrature rule

                        local u_exact = density_wave.(xq, yq, sol.t[end], equations)
                        local error = Vq * parent(u) - u_exact
                        local wJq = rd.wq.*rd.Vq*md.J
                        local num = sqrt(sum(wJq .* norm.(error).^2))
                        local den = sqrt(sum(wJq .* norm.(u_exact).^2))
                        local L2_rel_norm = num / den
                        local L2_norm = num

                        local folder = "Density_wave_Data/"
                        # local filename = "N$(N)_K$(K1D)_$(DG_type)_NO_AV.jld2"
                        local filename = "N$(N)_K$(K1D)_$(DG_type)_$(AV_discretization).jld2"
                        @save joinpath(folder, filename) md rd sol L2_norm L2_rel_norm


                end
        end

end

#shu_isentropic_euler_vortex Problem
if initial_condition_type == :shu_isentropic_vortex

        N_list = [3,4]
        K_list = [128]

        DG_type = :modal
        AV_discretization = :LDG

        for N in N_list
                for K1D in K_list
                        #Mesh Parameters
                        # N = 1 #Order of Polynomial
                        # K1D = 8 #Number of elements in each direction
                        display("N = $(N), K = $(K1D)")
                        if DG_type == :nodal
                                local rd = RefElemData(Tri(), SBP(), N)
                        elseif DG_type == :modal
                                local rd = RefElemData(Tri(), N)
                        end

                        local VXY, EToV = uniform_mesh(rd.element_type, K1D)
                        local L = 10.0 #Length of the domain
                        local VXY = VXY.*L #Scale the mesh
                        local md = MeshData((VXY,EToV), rd; is_periodic = true)
                        local x,y = md.x, md.y
                        local u = rd.Pq * shu_isentropic_vortex.(md.xq,md.yq,0.0,L,gamma,equations)
                        local Final_time = 15.0
                        #local Final_time = sqrt(2)*5.0/sqrt(2/gamma)
                        #local Final_time = 2.0*sqrt(gamma)*L - 5.0
                        local tspan = (0.0, Final_time)
                        local du = similar(VectorOfArray(u))
                        local params = (; rd, md, equations, interface_flux = flux_hllc, AV_discretization, gamma)
                        local ode = ODEProblem(rhs!, VectorOfArray(u), tspan, params)
                        local sol = solve(ode, Tsit5(); abstol = 1e-10, reltol = 1e-8, 
                                        saveat=LinRange(tspan..., 50), 
                                                callback=AliveCallback(alive_interval=10))
                        local u = parent(sol.u[end])


                        # higher accuracy quadrature rule
                        local rq, sq, wq = quad_nodes(rd.element_type, rd.N + 5)
                        local Vq = vandermonde(rd.element_type, rd.N, rq, sq) / rd.VDM
                        local xq, yq = Vq * md.x, Vq * md.y

                        # local xp, yp, up = rd.Vp * x, rd.Vp * y, rd.Vp * getindex.(u,1)

                        # display(scatter(vec(md.x), vec(md.y), zcolor=vec(getindex.(u,1)), msw=0, ms = 1, legend=false, ratio=1, cam=(0,90),
                        #         title="DG_Density", xlabel="x", ylabel="y", c=:red))

                        local alpha = pi/4
                        local M_inf = sqrt(2.0/gamma)
                        local sigma = 1.0
                        local R = 1.0
                        local Beta = M_inf*5.0*sqrt(2.0)*exp(0.5)/(4*pi)
                        local ux = M_inf * cos(alpha)
                        local vx = M_inf * sin(alpha)

                        local u_exact = shu_isentropic_vortex.(xq .- ux*sol.t[end], yq .- vx*sol.t[end], 0.0, L, gamma, equations)
                        # local u_exact = shu_isentropic_vortex.(xq, yq, sol.t[end], L, gamma, equations)
                        local error = Vq * u - u_exact
                        # up_exact = rd.Vp * getindex.(u_exact,1)

                        # display(scatter(vec(xp), vec(yp), zcolor=vec(up_exact), msw=0, ms = 1, legend=false, ratio=1, cam=(0,90),
                        # title="exact Density", xlabel="x", ylabel="y", c=:red))

                        # L2_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * u - u_exact)).^2)
                        local wJq = wq.*Vq*md.J
                        local L2_norm = sqrt(sum(wJq .* norm.(error).^2))

                        local L2_rel_norm = L2_norm / sqrt(sum(wJq .* norm.(u_exact).^2))

                        display("N = $(N), L2 norm = $(L2_norm)")

                        local base_path = @__DIR__

                        local folder = joinpath(base_path, "Vortex_Data_10")
                        mkpath(folder)
                        local filename = "N$(N)_K$(K1D)_T$(sol.t[end])_$(DG_type)_$(AV_discretization)_$(sol.prob.p.interface_flux).jld2"
                        
                        @save joinpath(folder, filename) md rd sol u_exact L2_norm L2_rel_norm
                        
                end
        end
end

# Riemann problem Configuration 1
if initial_condition_type == :Riemann_problem_config_1
        VXY, EToV = uniform_mesh(rd.element_type, K1D)
        md = MeshData((VXY,EToV), rd; is_periodic = true)
        x,y = md.x, md.y
        u = rd.Pq * Riemann_problem_config_1.(md.xq,md.yq,equations)
        Final_time = 0.25
        tspan = (0.0, Final_time)
        du = similar(VectorOfArray(u))
        params = (; rd, md, equations, interface_flux = flux_lax_friedrichs, AV_discretization)
        ode = ODEProblem(rhs!, VectorOfArray(u), tspan, params)
        sol = solve(ode, SSPRK43(); abstol = 1e-8, reltol = 1e-6, 
        saveat=LinRange(tspan..., 50), 
        callback=AliveCallback(alive_interval=10))
        u = parent(sol.u[end])
        xp, yp, up = rd.Vp * x, rd.Vp * y, rd.Vp * getindex.(u,1)
        display(scatter(vec(xp), vec(yp), zcolor=vec(up), msw=0, ms = 0.5, legend=false, ratio=1, cam=(0,90),
                title="DG_Density", xlabel="x", ylabel="y", c=:red))
        if DG_type == :nodal && AV_discretization == :LDG
                @save "Riemann_config_1_data/N$(N)_K$(K1D)_T$(Final_time)_nodal_LDG.jld2" md rd sol 
        elseif DG_type == :nodal && AV_discretization == :BR1
                @save "Riemann_config_1_data/N$(N)_K$(K1D)_T$(Final_time)_nodal_BR1.jld2" md rd sol 
        elseif DG_type == :modal && AV_discretization == :BR1
                @save "Riemann_config_1_data/N$(N)_K$(K1D)_T$(Final_time)_modal_BR1.jld2" md rd sol 
        elseif DG_type == :modal && AV_discretization == :LDG           
                @save "Riemann_config_1_data/N$(N)_K$(K1D)_T$(Final_time)_modal_LDG.jld2" md rd sol 
        end
end

