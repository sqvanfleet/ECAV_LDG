using StartUpDG
using LinearAlgebra 
using StaticArrays
using OrdinaryDiffEqSSPRK, RecursiveArrayTools
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

        if all(md.is_periodic .== (false,false)) 

                # --- Boundary Conditions Setup ---
                tol = 100 * eps()
                
                # Identify boundary face indices based on physical coordinates
                xb = md.xf[md.mapB]
                yb = md.yf[md.mapB]

                top_bottom_boundary = findall(@. abs(yb - 1.0) < tol || abs(yb) < tol)
                left_boundary       = findall(@. abs(xb) < tol)
                # right_boundary      = findall(@. abs(xb - 2.0) < tol)

                # 1. Top and Bottom Reflective Slip-Wall BCs
                for idx in top_bottom_boundary
                        b_node = md.mapB[idx]
                        rho, v1, v2, p = cons2prim(uM[b_node], equations)
                        nx, ny = md.nx[b_node], md.ny[b_node]
                        
                        # Normal velocity component: v_n = v1*nx + v2*ny
                        v_normal = v1 * nx + v2 * ny
                        
                        # Reflect normal velocity component
                        v1_reflect = v1 - 2 * v_normal * nx
                        v2_reflect = v2 - 2 * v_normal * ny
                        
                        uP[b_node] = prim2cons(SVector(rho, v1_reflect, v2_reflect, p), equations)
                end

                # 2. Inflow BC (Left Boundary, x = 0)
                # Prescribe exact incoming left-state/vortex flow
                for idx in left_boundary
                        b_node = md.mapB[idx]
                        uP[b_node] = Shock_vortex_interaction(md.xf[b_node], md.yf[b_node], equations, gamma)
                end

                # # 3. Outflow BC (Right Boundary, x = 2)
                # # Fix the unperturbed downstream post-shock state
                # for idx in right_boundary
                #         b_node = md.mapB[idx]
                #         uP[b_node] = Shock_vortex_interaction(md.xf[b_node], md.yf[b_node], equations, gamma)
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