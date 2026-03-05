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

data = load("Shock_vortex_interaction_Data/N3_K_64x32_T0.7_modal_LDG.jld2", "md", "rd", "sol")

rd = data[1]
md = data[2]
sol = data[3]

function calc_grad(u, equations, rd, md)
    drhodx = (md.rxJ .* (rd.Dr * rho) + md.sxJ .* (rd.Ds * rho)) ./ md.J
    drhody = (md.ryJ .* (rd.Dr * rho) + md.syJ .* (rd.Ds * rho)) ./ md.J
    g = sqrt.((rd.Vp * drhodx).^2 + (rd.Vp * drhody).^2)
end
    

