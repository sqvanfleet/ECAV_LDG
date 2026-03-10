# On the choice of viscous discontinuous Galerkin discretization for entropy correction artificial viscosity methods
This repository includes the Julia code for each numerical example given in Section 5 of 
[this paper](https://arxiv.org/abs/2602.23210).  The goal is that any reader may reproduce the results shown for each example.  

## Collaborator 
- Jesse Chan (https://www.ae.utexas.edu/people/faculty/faculty-directory/chan)

## Contents 

There are three directories in this repository, seven of them are aggregation-diffusion equation examples, which correspond to the examples in Section 5:
1. `Burgers_2D` (example 5.1)
2. `EulerEq_1D` (example 5.4, 5.6, and 5.7)
3. `EulerEq_2D` (example 5.3, 5.5) 

## Example 5.1 2D Inviscid Burgers Equation, Figure 1
1. Run `Burgers_2D.jl`.
2. Run `make_plots.jl`, plots will be generated in the `Plots` folder.

## Example 5.3 Shu Isentropic Euler Vortex, Table 1
1. run `Shu_Isentropic_Vortex.jl` to generate the LaTeX code for Table 1 in the output.

## Example 5.4 Stationary Contact Wave and Contact Preservation, Figure 2
1. run `stationary_contact_wave_piecewise.jl` to generate Figure 2a in the `Contact_wave_data` folder.
2. run `stationary_contact_wave_piecewise_smooth.jl` to generate Figure 2b in the `Contact_wave_data` folder.

## Example 5.5 Shock-vortex interaction, Figure 3
1. run the `shock_vortex_interaction.jl` to generate the `.png` files for Figure 3 in the `Shock_vortex_interaction_plots` folder.
2. the `.vtu` files can be used in Paraview to produce Figure 3.

## Example 5.6 1D Density Wave Figure 4, Figure 5, Figure 6, and Table 2
1. 

## Example 5.7 Shu-Osher problem Figure 7 and Figure 8
1. Navigate to the `EulerEq_1D` folder
2. open `1D Euler AV.jl` and uncomment line 306 and comment out lines 302 through 305 and comment out line 307
3. on line 291 change `N = 3` and `K1D = 100` and run `EulerEq_1D`
4. the plots in Figure 7 and Figure 8 (a) are in the `Shu_Osher_plots` folder
5. in `1D Euler AV.jl` comment out line 285 and uncomment line 286
6. the plot in Figure 8 (b) is in the `Shu_Osher_plots` folder
