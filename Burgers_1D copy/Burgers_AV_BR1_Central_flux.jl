using StartUpDG
using OrdinaryDiffEq
using Sundials
using JLD2

N = 3
K1D = 20 
rd = RefElemData(Line(),SBP(), N)
md = MeshData(uniform_mesh(Line(), K1D), rd;
              is_periodic=true)
f(u) = u^2/2
psi(u) = u^3/6

function rhs!(du, u, params, t)
    (; rd, md) = params
    invMQTr = -rd.M \ (rd.Dr' * rd.M) #M^(-1)*Q^(T) Why is Q^(T) = DG derivative times the mass matrix?
    uM = rd.Vf * u #Face quadrature points of u, in otherwords, this is u at the boundary on the interior of the element
    uP = uM[md.mapP] #Gets the values of u on the neighboring elements.
    lambda = @. max(abs(uP), abs(uM)) # gather the largest abs of f'(u) = u in this case for the upwind flux
    # interface_flux = @. 0.5 * (f(uP) + f(uM)) * md.nxJ #- 0.5 * lambda * (uP - uM)  #Upwind flux
    interface_flux = @. 1/6 * (uP^2 + uM * uP + uM^2) * md.nxJ - 0.0 * lambda * (uP - uM) 

    du .= md.rxJ .* (invMQTr * f.(u)) + rd.LIFT * interface_flux # Right hand side for the, LIFT matrix has Minv built into it.

    #BR1
    beta = 0
    #LDG
    #beta = @. sign(md.nxJ)

    #sigma = (md.rxJ .* (invMQTr * u) + rd.LIFT * (@. 0.5 * (uP + uM))) ./ md.J
    sigma = (md.rxJ .* (invMQTr * u) + rd.LIFT * (@. 0.5 * (uP + uM) * md.nxJ - 0.5 * beta * (uP - uM) * md.nxJ)) ./ md.J

    K1D = md.num_elements
    delta = zeros(K1D)
    num = zeros(K1D)
    den = zeros(K1D)
    epsilon = zeros(K1D)
    dvdxJ = md.rxJ .* (rd.Dr * u)
    term1 = @. -rd.wq * f(u) * dvdxJ
    term2 = @. rd.wf * psi(uM) * md.nxJ
    #den1 = sum(rd.M * sigma.^2 .* md.J,dims=1)
    for i in 1:K1D
        delta[i] = sum(term1[:,i]) + sum(term2[:,i])
        # Comput the numerator of eq (22)
        num[i] = -min(0, delta[i])
        # Compute the denomenator of eq (22)
        den[i] = sum(rd.wq .* sigma[:,i].^2 .* md.J[:,i])
        epsilon[i] = num[i] * den[i] / (1e-14 + den[i]^2)
        sigma[:,i] .*= epsilon[i] 
    end

    #@show extrema(epsilon), t

    sigmaM = rd.Vf * sigma
    sigmaP = sigmaM[md.mapP]
    du .-= md.rxJ .* (invMQTr * sigma) + rd.LIFT * (@. 0.5 * (sigmaP + sigmaM) * md.nxJ + 0.5 * beta * (sigmaP - sigmaM) * md.nxJ)
    @. du /= -md.J

    return du
end
x = md.x
u = @. 0.5 - sin(pi * x)
params = (; rd, md)
tspan = (0, 0.365)
ode = ODEProblem(rhs!, u, tspan, params)
sol = solve(ode, SSPRK43(), abstol = 1e-6, reltol = 1e-4, 
            saveat=LinRange(tspan..., 50))


subdirectory = "Data"
filename = "Central_flux_BR1"

# Check if the subdirectory exists, and create it if necessary
if !isdir(subdirectory)
    mkpath(subdirectory)  # Create the subdirectory if it doesn't exist
end

@save joinpath(subdirectory, filename) sol

#Compute epsilon using the solution at each time interval.

sol_size = size(sol)
epsilon_evolution = zeros(K1D,sol_size[3])
delta = zeros(K1D)
num = zeros(K1D)
den = zeros(K1D)


for (j,u) in enumerate(sol.u)
    invMQTr = -rd.M \ (rd.Dr' * rd.M)

    #BR1
    beta = 0
    #LDG
    #beta = @. sign(md.nxJ)
    
    uM = rd.Vf * u
    uP = uM[md.mapP]
    #LDG
    sigma = (md.rxJ .* (invMQTr * u) + rd.LIFT * (@. 0.5 * (uP + uM) * md.nxJ - 0.5 * beta * (uP - uM) * md.nxJ)) ./ md.J
    
    
    dvdxJ = md.rxJ .* (rd.Dr * u)
    term1 = @. -rd.wq * f(u) * dvdxJ
    term2 = @. rd.wf * psi(uM) * md.nxJ
    for i in 1:K1D
        delta[i] = sum(term1[:,i]) + sum(term2[:,i])
        # Comput the numerator of eq (22)
        num[i] = -min(0, delta[i])
        # Compute the denomenator of eq (22)
        den[i] = sum(rd.wq .* sigma[:,i].^2 .* md.J[:,i])
        epsilon_evolution[i,j] = num[i] * den[i] / (1e-14 + den[i]^2)
    end

end


#=
anim = @animate for i in 1:length(sol.u)
    u = sol.u[i]
    #.Vp adds more nodes for plotting
    str1 = string(sol.t[i])
    str2 = "BR1_solution at time" * " " * str1
    scatter(rd.Vp * x, rd.Vp * u, legend=false, ylims=(-1, 2), title=str2)
end 

gif(anim,"Gifs/AV_BR1_Central_flux_solution_evolution.gif",fps=10)

#=
@gif for i in 1:sol_size[3]
    time = sol.t[i]
    epsilon_gif = copy(epsilon_evolution[:,i]);
    scatter(x[1,:],epsilon_gif, legend=false,
        title ="epsilon evolution at time "*string(time), ylims=(0.0,maximum(epsilon_evolution)))
end filename="my_gif.gif" fps=30 #Save every frame, and set fps to 30
=#

plot(sol.t,transpose(maximum(epsilon_evolution,dims=1)),
    label = "", title="Burgers 1D AV BR1 Central Flux", xlabel="Time", ylabel="Maximum Epsilon Value")

# Save the plot to a file
pathname = "/Users/samvanfleet/Documents/Rice Artificial Viscosity/Burgers_1D/\
    Plots/Burgers_1D_AV_BR1_Central_flux_Epsilon_Evolution.pdf"
savefig(pathname)  # You can change the file format to .pdf, .svg, etc.
=#