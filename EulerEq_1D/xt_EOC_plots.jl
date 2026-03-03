using JLD2
using StartUpDG
using LinearAlgebra 
using StaticArrays
using Plots
using Printf
using Trixi
using Trixi.ForwardDiff
using OrdinaryDiffEq, RecursiveArrayTools
using LaTeXStrings


psi(u, normal, ::CompressibleEulerEquations1D) = u[2] * normal
dudv(v, equations) = ForwardDiff.jacobian(v -> entropy2cons(v, equations), v)


function initial_condition_density_wave(x, equations::CompressibleEulerEquations1D; 
                                        amplitude = 0.5)
    rho = 1.0 + amplitude * (sin(2*pi*x)+(abs(x) < 0.5))
    #rho = 1.0 + amplitude * (sin(2*pi*x))
    u = 0.0
    p = 1.0

    return prim2cons(SVector(rho, u, p), equations)
end


gamma = 1.4
equations = CompressibleEulerEquations1D(gamma)


#Load the course mesh solution
data1 = load("Density_wave_Data/N3_K80.jld2", "md", "rd", "sol")

md = data1[1]
rd = data1[2]
sol = data1[3]

params = (; rd, md, equations, interface_flux = flux_hllc)

x = md.x

xbound = [md.xf[1,:],md.xf[2,end]] # make a list of element boundaries

xcenter = md.xf[1,:] .+ (md.xf[end,end] - md.xf[1,1])/(2*md.K) # make a list of element centers

xt_error = zeros(length(sol.t), size(xcenter,1))
xt_rel_error = zeros(length(sol.t), size(xcenter,1))

for (i,u) in enumerate(sol.u)
    t = sol.t[i]
    local u_exact = initial_condition_density_wave.(md.xq , equations)
    error = sqrt.(sum(md.wJq.* norm.(rd.Vq*parent(u)-u_exact).^2, dims = 1))
    den = sqrt.(sum(md.wJq.* norm.(u_exact).^2, dims = 1))
    rel_error = error./den
    xt_error[i,:] = error
    xt_rel_error[i,:] = rel_error
end

#compute L2 error over the entire domain at the final time
u_exact_course = initial_condition_density_wave.(md.xq .- 0.1 * sol.t[end], equations)
L2_error_course = sqrt(sum(md.wJq .* norm.(rd.Vq * parent(sol.u[end]) - u_exact_course).^2))

#load the dense mesh solution
data2 = load("Density_wave_Data/N3_K160.jld2", "md", "rd", "sol")


md = data2[1]
rd = data2[2]
sol = data2[3]

xt_error_dense = zeros(length(sol.t), size(xcenter,1))
xt_rel_error_dense = zeros(length(sol.t), size(xcenter,1))

for (i,u) in enumerate(sol.u) 
    t = sol.t[i]
    local u_exact = initial_condition_density_wave.(md.xq .- 0.1 * t, equations)
    error = sqrt.(sum(md.wJq.* norm.(rd.Vq*(parent(u))-u_exact).^2, dims = 1)) #Fix this error sqrt + sqrt 
    error_1 = error[:,1:2:end]
    error_2 = error[:,2:2:end]
    error_dense = error_1 + error_2
    den = sqrt.(sum(md.wJq.* norm.(u_exact).^2, dims = 1))
    den_1 = den[:,1:2:end]
    den_2 = den[:,2:2:end]
    den_dense = den_1 + den_2
    rel_error = error_dense ./ den_dense
    xt_error_dense[i,:] = error_dense
    xt_rel_error_dense[i,:] = rel_error
end

EOC = log.(xt_error./xt_error_dense) ./ log(2)


EOC_plot = heatmap(xcenter, sol.t, EOC,camera = (0,90),colorbar = true,
        xlabel = L"$x$", ylabel = L"$t$", title = "EOC for discontinuous initial data",
        rightmargin = 10Plots.mm)

Error_plot = heatmap(xcenter, sol.t, xt_error, camera = (0,90),colorbar = true,
        xlabel = L"$x$", ylabel = L"$t$", title = "Error (Semi-Log y) for discontinuous initial data",
        rightmargin = 10Plots.mm)



display(EOC_plot)

display(Error_plot)

#compute error on entire domain at final time

u_exact_dense = initial_condition_density_wave.(md.xq .- 0.1 * sol.t[end], equations)
L2_error_dense = sqrt(sum(md.wJq .* norm.(rd.Vq * parent(sol.u[end]) - u_exact_dense).^2))

error_ratio = L2_error_course/L2_error_dense
EOC_final_time = log.(error_ratio) ./ log(2)


display(EOC_final_time)
display(L2_error_dense)
display(L2_error_course)





