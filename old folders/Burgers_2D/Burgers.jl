using StartUpDG
using OrdinaryDiffEq
using Plots
using Trixi: AliveCallback
using JLD2

N = 3
K1D = 20
rd = RefElemData(Tri(), SBP(), N)
md = MeshData(uniform_mesh(Tri(), K1D), rd; is_periodic=true)

f1(u) = u^2/2
f2(u) = u^2/2

psi1(u) = u^3/6
psi2(u) = u^3/6


function rhs!(du, u, params, t)
    (; rd, md) = params
    invMQTr = -rd.M \ (rd.Dr' * rd.M)
    invMQTs = -rd.M \ (rd.Ds' * rd.M)
    uM = rd.Vf * u
    uP = uM[md.mapP]
    interface_flux = @. 1/2 * (f1(uP) + f1(uM)) * md.nxJ + 
                        1/2 * (f2(uP) + f2(uM)) * md.nyJ - 0.5 * (uP - uM)

    du .= md.rxJ .* (invMQTr * f1.(u)) + md.sxJ .* (invMQTs * f1.(u)) + 
          md.ryJ .* (invMQTr * f2.(u)) + md.syJ .* (invMQTs * f2.(u)) + rd.LIFT * interface_flux

    
    beta =  @. sign(md.nxJ + md.nyJ)

    #beta = 0

    sigma_1 = ((md.rxJ + md.ryJ) .* (invMQTr * u) + rd.LIFT * (@. 0.5*(uP + uM)*md.nxJ - 0.5 * beta * (uP - uM) * md.nxJ)) ./ md.J
    sigma_2 = ((md.sxJ + md.syJ) .* (invMQTs * u) + rd.LIFT * (@. 0.5*(uP + uM)*md.nyJ - 0.5 * beta * (uP - uM) * md.nyJ)) ./ md.J
    
    #epsilon = 0.001

    K = md.num_elements
    
    delta = zeros(K)
    num = zeros(K)
    den = zeros(K)
    epsilon = zeros(K)

    #dvdxJ = -(md.rxJ + md.ryJ) .* (rd.Dr' * rd.M * u)
    #dvdyJ = -(md.sxJ + md.syJ) .* (rd.Ds' * rd.M * u)
    dvdxJ = -md.rxJ .* (rd.Dr' * rd.M * u) - md.sxJ .* (rd.Ds' * rd.M * u) 
    dvdyJ = -md.ryJ .* (rd.Dr' * rd.M * u) - md.syJ .* (rd.Ds' * rd.M * u)
    term1 = @. f1(u)*dvdxJ + f2(u)*dvdyJ
    term2 = @. rd.wf * psi1(uM) * md.nxJ + rd.wf * psi2(uM * md.nyJ)

    delta = sum(term1, dims=1) + sum(term2, dims=1)
    num = @. -min(0, delta)
    den = sum(rd.M * (sigma_1.^2 + sigma_2.^2) .* md.J, dims=1)
    epsilon = @. num*den/(1e-14 + den^2)

    
    sigma_1 = sigma_1.*epsilon
    sigma_2 = sigma_2.*epsilon

    
    sigma_1M = rd.Vf * sigma_1
    sigma_1P = sigma_1M[md.mapP]
    
    sigma_2M = rd.Vf * sigma_2
    sigma_2P = sigma_2M[md.mapP]
    
    du .-= (md.rxJ + md.ryJ) .* (invMQTr * sigma_1) + (md.sxJ + md.syJ) .* (invMQTs * sigma_2) +
            rd.LIFT * (@. 0.5 * (sigma_1P + sigma_1M) * md.nxJ - 0.5 * beta * (sigma_1P - sigma_1M) * md.nyJ) +
            rd.LIFT * (@. 0.5 * (sigma_2P + sigma_2M) * md.nyJ - 0.5 * beta * (sigma_2P - sigma_2M) * md.nyJ)  

    @. du /= -md.J

end

function Gaussian_initial_condition(x,y)
    return exp((-(1.0*x)^2 - (1.0*y)^2))
end

function piecewise_initial_condition(x,y)
    if x^2 + y^2 < 0.5
        return 0.35*pi
    else
        return 0.25*pi
    end
end

x = md.x
y = md.y
u = @. Gaussian_initial_condition(x,y)

params = (; rd, md)
tspan = (0, 1.0)
ode = ODEProblem(rhs!, u, tspan, params)
sol = solve(ode, SSPRK43(); abstol = 1e-6, reltol = 1e-4, saveat=LinRange(tspan..., 50), 
            callback=AliveCallback(alive_interval = 100))

@save "Data/solution.jld2" sol

u = sol.u[end]
xp, yp, up = rd.Vp * x, rd.Vp * y, rd.Vp * u
scatter(vec(xp), vec(yp), zcolor=vec(up), msw=0, legend=false, 
        ms=1, ratio=1, cam=(0,90), title="Final Solution")

#=            
str1 = "Advection_2d_no_AV"

xp, yp, up = rd.Vp * x, rd.Vp * y, rd.Vp * u

anim = @animate for i in eachindex(sol.u)
    u = sol.u[i]
    up = rd.Vp * u
    #.Vp adds more nodes for plotting
    str2 = string(sol.t[i])
    str3 = str1 * " " * str2
#    scatter(xp, yp, up, zcolor=up, msw=0, legend=false, ratio=1, cam=(0,90),title = str3)
    scatter(vec(xp), vec(yp), zcolor=vec(up), msw=0, legend=false, ratio=1, title = str3)
end 

pathname = str1 * ".gif"

gif(anim,pathname,fps=10)




#surface(x,y,u)

xp, yp, up = rd.Vp * x, rd.Vp * y, rd.Vp * u
#scatter(xp, yp, up, zcolor=up, msw=0, legend=false, ratio=1, cam=(0,90))
=#