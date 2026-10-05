# paper-artificial-viscosity-entropy-stable-2025

This repository contains information and code to reproduce the results presented in the preprint ["An artificial viscosity approach to high order entropy stable discontinuous Galerkin methods"](https://arxiv.org/abs/2501.16529). 

You will need to first activate and instantiate the Julia package environment to run any code. These codes were tested on Julia 1.10.7. 

# Paper codes

The codes are listed in the order the corresponding figures appear in the preprint. 

* `spurious_nullspaces.jl`: computes and plots spurious null space modes of the 1D and 2D BR-1 gradients.
* `density_wave_2D_convergence.jl`: computes L2 errors and convergence rates for the 2D density wave. 
* `density_wave_over_time_1D.jl`: computes L2 error over time for the 1D density wave
* `linear_stability_1D.jl`: computes the Jacobian and spectra using ForwardDiff.jl of a linearization around the 1D density wave
* `entropyresidual_viscosity_2D_convergence.jl`: computes the maximum entropy residual and viscosity coefficient over a sequence of refined meshes. 
* `AV_1D.jl`: solver for the 1D modified Sod and Shu-Osher problems using artificial viscosity-based entropy stable DG.
* `flux_diff_1D.jl`: solver for the 1D modified Sod and Shu-Osher problems using flux differencing entropy stable DG formulations.
* `AV_2D.jl`: solver for the 2D Riemann problem and long-time Kelvin-Helmholtz instability using artificial viscosity-based entropy stable DG.

## Appendix codes
* `compare_1D_elementwise_subcell_AV.jl`: compares element-wise and subcell artificial viscosity coefficients for the Shu-Osher problem. 
* `entropy_correction_2D.jl`: implements a modified version of the local entropy correction term from Abgrall, Offner, and Ranocha 2018.
