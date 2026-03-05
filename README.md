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
1. Navigate to the `Burgers_2D` folder
2. Run `Burgers_2D.jl`
3. Run `make_plots.jl`, the plots in Figure 1 will be in the Plots folder

## Example 5.3 Shu Isentropic Euler Vortex, Table 1
1. navigate to the `EulerEq_2D` folder
2. open the `2D_Euler_AV.jl` file and ensure line 142 is
   uncommented and lines 121,143,144 are commented, then run `2D_Euler_AV.jl`
3. run `Vortex_Problem_Convergence_Tabel.jl that outputs the LaTeX code for Table 1.

## Example 5.4 Stationary Contact Wave and Contact Preservation, Figure 2
1. navigate to the `EulerEq_1D` folder
2. open the `1D Euler AV.jl` file and ensure line 305 is uncommented and lines 302,303,304, 306,307 are commented.
3. on line 292 ensure `K1D = 80` and run `1D Euler AV.jl` for `N = 2`, `N = 4`, `N=6`, and `N = 8` on line 291
4. run `Contact_wave_error_evolution.jl` to generate Figure 2 (b).
5. open `initial_conditions.jl` and comment line 23 and uncomment lines 24 through 28.
6. repeat steps 1 through 4 to generate Figure 2 (a)

## Example 5.5 Shock-vortex interaction, Figure 3
1. navigate to the `Euler_2D` folder
2. open the `2D_Euler_AV.jl` file and ensure line 144 is not commented out and that lines 14 though 143 are commented out.
3. on line 151 change `N = 2` and on line 152 change `K1D = 64` and run the `2D_Euler_AV.jl` file.
4. The .vtu files in the `Shock_vortex_interaction_plots` can be used by Paraview to creat the plots in Figure 3

## Example 5.6 1D Density Wave Figure 4, Figure 5, Figure 6, and Table 2
1. Navigate to the `EulerEq_1D` folder
2. open `1D Euler AV.jl` and uncomment line 302 and comment out lines 303 through 307
3. on line 291 change `N = 5` and on line 292 change `K1D = 8` and run `1D Euler AV.jl`,
   the first row of Table 2 is in the output
5. on line 292 change `K1D = 16` and run `1D Euler AV.jl` and the second row of Table 2 is in the output
7. The plots in Figure 4, 5, and 6 are now in the `Density_Wave_Plots` folder
8. on line 292 change `K1D = 24` and run `1D Euler AV.jl` and the third row of Table 2 is in the output

## Example 5.7 Shu-Osher problem Figure 7 and Figure 8
1. Navigate to the `EulerEq_1D` folder
2. open `1D Euler AV.jl` and uncomment line 306 and comment out lines 302 through 305 and comment out line 307
3. on line 291 change `N = 3` and `K1D = 100` and run `EulerEq_1D`
4. the plots in Figure 7 and Figure 8 (a) are in the `Shu_Osher_plots` folder
5. in `1D Euler AV.jl` comment out line 285 and uncomment line 286
6. the plot in Figure 8 (b) is in the `Shu_Osher_plots` folder
