using StartUpDG
using Plots

rd = RefElemData(Line(), N=2)
md = MeshData(uniform_mesh(Line(), 8), rd; is_periodic=true)

function rhs(u, rd, md)
    uM = rd.Vf * u
    uP = uM[md.mapP]
    sigma = md.rxJ .* (rd.Dr * u) + rd.LIFT * (@. 0.5 * (uP - uM) * md.nxJ)
    sigmaM = rd.Vf * sigma
    sigmaP = sigmaM[md.mapP]
    return md.rxJ .* (rd.Dr * sigma) + rd.LIFT * (@. 0.5 * (sigmaP - sigmaM) * md.nxJ)
end

A = Matrix(LinearMap(u -> rhs(reshape(u, :, md.num_elements), rd, md), length(md.x)))
M = kron(Diagonal(md.J[1,:]), rd.M)

K = M * A
@show norm(K - K') < 100 * eps()
lam, W = eigen(Symmetric(K))

Vp = vandermonde(Line(), rd.N, LinRange(-1,1,100)) / rd.VDM
plt = plot(pad_nans(Vp * md.x), pad_nans(Vp * reshape(W[:,end-1], :, md.num_elements)), 
           linewidth=2, leg=false)#, xformatter=:none, yformatter=:none)
png(plt, "spurious_modes_1d.png")

# 2D null spaces
rd = RefElemData(Tri(), N=2)
md = MeshData(uniform_mesh(Tri(), 4), rd; is_periodic=true)

function rhs(u, rd, md)
    uM = rd.Vf * u
    uP = uM[md.mapP]
    sigma_x = md.rxJ .* (rd.Dr * u) + md.sxJ .* (rd.Ds * u) + 
        rd.LIFT * (@. 0.5 * (uP - uM) * md.nxJ)
    sigma_y = md.ryJ .* (rd.Dr * u) + md.syJ .* (rd.Ds * u) + 
        rd.LIFT * (@. 0.5 * (uP - uM) * md.nyJ)

    sigma_x_M, sigma_y_M = rd.Vf * sigma_x, rd.Vf * sigma_y
    sigma_x_P, sigma_y_P = sigma_x_M[md.mapP], sigma_y_M[md.mapP]
    d_sigma_x = md.rxJ .* (rd.Dr * sigma_x) + md.sxJ .* (rd.Ds * sigma_x) + 
        rd.LIFT * (@. 0.5 * (sigma_x_P - sigma_x_M) * md.nxJ)
    d_sigma_y = md.ryJ .* (rd.Dr * sigma_y) + md.syJ .* (rd.Ds * sigma_y) + 
        rd.LIFT * (@. 0.5 * (sigma_y_P - sigma_y_M) * md.nyJ)
    return d_sigma_x, d_sigma_y
end

A = Matrix(LinearMap(u -> vcat(vec.(rhs(reshape(u, :, md.num_elements), rd, md))...), 
           2 * length(md.x), length(md.x)))
M = kron(Diagonal(md.J[1,:]), rd.M)

K = kron(I(2), M) * A
_, s, W = svd(K)

plt = scatter(pad_nans(Vp * md.x), pad_nans(Vp * md.y), 
              zcolor=pad_nans(Vp * reshape(W[:,end-1], :, md.num_elements)), 
              msw=0, ms=1, leg=false, ratio=1, axis=([], false),
              xformatter=:none, yformatter=:none)
png(plt, "spurious_modes_2d.png")
