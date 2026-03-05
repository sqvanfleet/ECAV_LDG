# On the choice of viscous discontinuous Galerkin discretization for entropy correction artificial viscosity methods
This repository includes the Julia code for each numerical example given in section 5 of 
[this paper](https://arxiv.org/abs/2602.23210).  The goal is that any reader may reproduce the results shown for each example.  

## Collaborator 
- Jesse Chan (https://www.ae.utexas.edu/people/faculty/faculty-directory/chan)

## Contents 

There are three directories in this repository, seven of them are aggregation-diffusion equation examples, which correspond to examples 4.1 - 4.4:
1. `Burgers_2D` (example 5.1)
2. `EulerEq_1D` (example 5.4, 5.6, and 5.7)
3. `EulerEq_2D` (example 5.3, 5.5)

## MATLAB requirements

All of the files in this repository were tested on MATLAB version R2022b.  The examples corresponding to the Landau equation require the [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html), but can be modified to be used in serial. 

## Project description

### `particle_method.m` or `particle_method_2d_parallel.m` scripts

#### Parameters, initialization, and mesh generation

- The `particle_method.m` or `particle_method_2d_parallel.m` scripts begin by defining several parameters including:
- `n_list` which is a matrix whose elements are how many particles you want to calculate the particle method solution with.  The entries of this array are chosen by the user.  The entries of `n_list` correspond to different values of $M$ from the paper, making the total number of particles $M^d$, where $d$ is the number of spatial dimensions.  The outermost loop in `particle_method.m` or `particle_method_2d_parallel.m` uses the loop variable `alpha` and loops over all values in `n_list`.  For each iteration of this loop the number of particles is
```matlab
n = n_list(alpha);
```
- For the 2D examples,
```matlab
Np = n^2
```
is the number of particles.  

- The porous medium example uses an additional parameter `m` which is the constant from the porous medium equation.

- Landau equation examples have two additional parameters `gamma` and `C_gamma`, which are parameters representing the collision kernel 
$$a_{ij}(\boldsymbol{x})=C|\boldsymbol{x}|^{\gamma}(|\boldsymbol{x}|^2\delta_{ij}-x_ix_j)$$

- `Xmax` or `Vmax` gives the computational domain which is the hypercube centered at the origin               
$[\mbox{Xmax},\mbox{Xmax}]^d$ or $[\mbox{Vmax},\mbox{Vmax}]^d$. These variables are chosen by the user.
- `dx` or `dv` is cell length so that the cell volume is $dx^d$ or $dv^d$, and is computed with
```matlab
dx = 2*Xmax/n
```
```matlab
dv = 2*Vmax/n
```
- The particle locations are initialized at the cell centers and for 1D examples we have
```matlab
x = (-Xmax+dx/2):dx:(Xmax-dx/2);
```
and for 2D examples we use
```matlab
v = (-Vmax+dv/2):dv:(Vmax-dv/2);
[vx,vy] = ndgrid(v);
Vx = vx(:);
Vy = vy(:);
```
to create two arrays `Vx` and `Vy` which contain the coordinates of the particles the $x$ and $y$ 
direction.
- `t0` is the initial time and chosen by the user.
- `f0` or `f_initial` are the initial conditions which are computed using the `exact.m`, `exact_2d.m`, or `non_bkw_initial_conditions.m` functions.  
- `w` are the particle weights and they are initialized using the midpoint rule.  In the 2D examples `W = w(:)` flattens the `w` array so that `W` is the same size as the `Vx` and `Vy` arrays. 
- To create the "blob solution" (equation 4.1), the reconstruction mesh is created.  For 1D examples
```matlab
Nr = n;
dxr = 2*Xmax/n; 
xr = (-Xmax+dxr/2):dxr:(Xmax-dxr/2); 
```
   
and for 2D examples
```matlab
Nr = n;
dvr = 2*Vmax/Nr;
vr = (-Vmax+dvr/2):dvr:(Vmax-dvr/2); 
[vrx,vry] = ndgrid(vr);
Vrx = vrx(:);
Vry = vry(:);
```
Creating this reference mesh is required to compute the errors or plot the blob solution.
- Once the mesh size `dx` or `dv` is choses the regularization parameter $\epsilon$ is chosen to be
```matlab
epsilon = 4*(0.4*(dv)^0.99)^2;
```
- The initial "reconstructed or blob solution" can then be calculated using the formula $$\sum^N_{p=1} w_p\varphi_{\varepsilon}(\boldsymbol{x} - \boldsymbol{x}_p).$$
####   The `psi_1d` and `psi_2d` functions
- The arguments for the `psi_1d` are
    - `x` represents the 1D mesh
    - `eps` which represents the regularization parameter
- The arguments for `psi_2d` are
    - `vx`, `vy`, represent the 2D mesh
    - `eps` represents the regularization parameter.
- These functions return an array the same size as `x` or `vx` and `vy`.  Each entry of the returned array is the mollifier function
  $$\varphi_{\varepsilon}(\boldsymbol{x}) = \frac{1}{2 \pi \varepsilon}\exp{\left(\frac{-|\boldsymbol{x}^2|}{2 \varepsilon}\right)},$$
  evaluated at the corresponding entry of `x` or `vx` and `vy` with the parameter `eps`.
  
```matlab
f = zeros(1,Nr);
for i = 1:Nr
    f(i) = sum(w.*psi_1d(xr(i)-x,epsilon));
end
```
and in 2D
```matlab
parfor i = 1:Nr
    for j = 1:Nr
        f(i,j) = sum(W.*psi_2d(vrx(i,j)-Vx,vry(i,j)-Vy,epsilon));
    end
end
```
- `tmax` is the final time and is chosen by the user.
- `dt` is the time step and is chosen by the user.
- `Nt` is the number of time steps
  ```matlab
  Nt = round((tmax-t0)/dt);
  ```
- `error_list` is an array that stores several quantities at each time step such as errors, mass, momentum, kinetic energy,
   energy, the dissipation term, fisher energy, and how many fixed point iterations were required for the relative error to
   be less than the given tolerance.
- `max_iter` is the maximum number of iterations the fixed point iteration can perform and is chosen by the user.
- `tol` is the tolerance used in the fixed point iteration and is chosen by the user.

#### For loop in time

- The for loop in time uses `nt` as the loop variable and goes from 1 to `Nt`
  
- `time` represents the current time
```matlab
time = t0+dt*nt
```

#### `right_hand_side` or `right_hand_side_parallel` functions

- The arguments for the `right_hand_side` function are
    - `w` represents the 1D wights 
    - `x` represents the 1D particle locations at the fixed point iteration `i`
    - `x_new` represents the particle locations at the fixed point iteration `i+1`
    - `xr` represents the 1D reference mesh
    - `dx` is the length of the reference mesh elements.  
    - `epsilon` is the regularization parameter
    - `n` is the number of particles
    - the `right_hand_side` function in the `particle_method_1D_porous_medium_discrete_gradient` uses an additional argument `m` and represents the constant from the porous medium equation.
      
- The output of `right_hand_side` is:
    -  array `gF` is the discrete gradient
  $$-\frac{1}{w_p}\overline{\nabla_{\boldsymbol{x}_p}E_A^{\varepsilon}}\left(\boldsymbol{X}^{n+1},\boldsymbol{X}^n\right) = \int_0^1 \nabla\_{\boldsymbol{x}_p} E\_{A}^{\varepsilon}\left(\boldsymbol{X}^n+s(\boldsymbol{X}^{n+1}-\boldsymbol{X}^n)\right)\mathrm{d}s,$$
  and is approximated with a 4 point Gauss-Legendre quadrature using the `lgwt` function.

- The arguments for the `right_hand_side_parallel` function are
    - `W` represents the particle weights
    - `Vx` represents the $x$ coordinates of the particle locations at the fixed point iteration `i`
    - `Vy` represents the $y$ coordinates of the particle locations at the fixed point iteration `i`
    - `Vx_new` represents the $x$ coordinates of the particle locations at the fixed point iteration `i+1`
    - `Vy_new` represents the $y$ coordinates of the particle locations at the fixed point iteration `i+1`
    - `Vrx` represents the $x$ coordinates of the 2D reference mesh 
    - `Vry` represents the $y$ coordinates of the 2D reference mesh 
    - `dv` represents the length of the reference mesh element making the element size `dv^2`
    - `epsilon` is the regularization parameter
    - `C_gamma` is the $C$ parameter from the Landau equation 
    - `gamma` is the $gamma$ parameter from the Landau equation
    - `Np` is the number of particles

- The outputs of `right_hand_side_parallel` are:
    - arrays `U_x` and `U_y`  are the discrete gradient
  $$-\frac{1}{w_p}\overline{\nabla_{\boldsymbol{x}_p}E_L^{\varepsilon}}\left(\boldsymbol{X}^{n+1},\boldsymbol{X}^n\right) = \int_0^1 \nabla\_{\boldsymbol{x}_p} E\_{L}^{\varepsilon}\left(\boldsymbol{X}^n+s(\boldsymbol{X}^{n+1}-\boldsymbol{X}^n)\right)\mathrm{d}s,$$
  and is approximated with a 4 point Gauss-Legendre quadrature using the `lgwt` function.

    - arrays `gF_x` and `gF_y` are the quantities used to calculate the fisher information
      $$F^{\varepsilon}(f^N) = \sum_{p = 1}^N\frac{1}{w_p}\left| \overline{\nabla_{\boldsymbol{x}_p} E_L^{\varepsilon}}(\boldsymbol{X}^{n+1},\boldsymbol{X}^n) \right|^2$$

    - `dissipation` is used to track the energy dissipation term
      
  $$D^{\varepsilon}(f^N) = -\frac{1}{2}\sum_{p,q=1}^N w_p w_q \left(\frac{\overline{\nabla_{\boldsymbol{x}\_p} E_L^{\varepsilon}}(\boldsymbol{X}^{n+1},\boldsymbol{X}^n)}{w_p} - \frac{\overline{\nabla_{\boldsymbol{x}\_q} E_L^{\varepsilon}}(\boldsymbol{X}^{n+1},\boldsymbol{X}^n)}{w_q}\right)A(\overline{\boldsymbol{x}\_p^n}-\overline{\boldsymbol{x}\_q^n})\left(\frac{\overline{\nabla_{\boldsymbol{x}\_p} E_L^{\varepsilon}}(\boldsymbol{X}^{n+1},\boldsymbol{X}^n)}{w_p} - \frac{\overline{\nabla_{\boldsymbol{x}\_q} E_L^{\varepsilon}}(\boldsymbol{X}^{n+1},\boldsymbol{X}^n)}{w_q}\right).$$


- The next for loop in the `particle_method` or `particle_method_2d_parallel` scripts is the loop that performs fixed point iteration that approximates the solution to the
  system resulting from the discrete gradient integrator.  This loop variable is `i` and goes from 1 to `max_iter` and represents the
  index of the fixed point iteration.

- Before this loop begins, the `right_hand_side` and `right_hand_side_parallel` functions are used to produce a
  forward Euler initial guess.  This works because of the discrete gradient property
  $$\overline{\nabla_{\boldsymbol{x}_p}E\_{A/L}^{\varepsilon}}\left(\boldsymbol{X}^{n},\boldsymbol{X}^n\right) = \nabla\_{\boldsymbol{x}_p}E^{\varepsilon}\_{A/L}(\boldsymbol{X}^n)$$

- For the 1D examples this is 
```matlab
dF = right_hand_side(w,x,x,xr,dx,epsilon,n);
xnew = x - dt*dF;
```
- For the 2D examples this is
```matlab
[U_x,U_y,gF_x,gF_y,dissipation] = right_hand_side_parallel(W,Vx,Vy,Vx,Vy,Vrx,Vry,dv,epsilon,C_gamma,gamma,Np);
Vx_new = Vx + dt*U_x;
Vy_new = Vy + dt*U_y;
```
- The fixed point iteration loop begins by defining `x_old` or `vx_old` and `vy_old` so that the relative error in the fixed point iteration can be calculated later on 
```matlab
x_old = x_new;
```
  or
```matlab
Vx_old = Vx_new;
Vy_old = Vy_old;
```

- Using the `right_hand_side` or `right_hand_side_parallel` functions, calculate overwrite `x_new` or `Vx_new` and `Vy_new` with a iteration of fixed point method
```matlab
gF = right_hand_side(w,x,xnew,xr,dx,epsilon,n);
xnew = x - dt*gF;
```
  or
```matlab
[U_x,U_y,gF_x,gF_y,dissipation] = right_hand_side_parallel(W,Vx,Vy,Vx_new,Vy_new,Vrx,Vry,dv,epsilon,C_gamma,gamma,Np);
Vx_new = Vx + dt*U_x;
Vy_new = Vy + dt*U_y;
```
    
- Using `x_old` and `x_new` (`Vx_old`, `Vy_old`, `Vx_new`, and `Vy_new` in the 2D exampels), the absolute error and relative error are computed and saved to `absres` and `relres` in 1D this is 
```matlab
absres = sqrt(sum((xnew - x_old).^2));
relres = absres/(sqrt(sum(xnew.^2)));
```
in 2D this is 

```matlab
absres = sqrt(sum((Vx_new - Vx_old).^2 + (Vy_new - Vy_old).^2));
relres = absres/(sqrt(sum(Vx_new.^2 + Vy_new.^2)));
```

- Using `relres` and `tol`, an `if` statement is used to determine if the relative error meets the stopping criterion.  If `relres` is less than `tol` the fixed point iteration loop breaks and the value of the index variable `i` is saved to the variable `numiter` which is later saved as an entry of the array `error_list`.  If the value of the index variable `i` ever is equal to `max_iter` then the time `time` and `relres` is saved in the array `Fixed_Point_list`
```matlab
if relres < tol
    disp(['fixed point iteration took ',num2str(i),' iterations at time ', num2str(time)])
    numiter = i;
    break
end
if i == max_iter
    Fixed_Point_list(nt,1) = time;
    Fixed_Point_list(nt,2) = relres;
    disp(['maximum number of iterations reached at time ',num2str(time)])
    numiter = i;
end
```

-  The fixed point iteration loop is complete and the particle locations are updated in 1D
```matlab
x = xnew;
```
and in 2D
```matlab
Vx = Vx_new;
Vy = Vy_new;
```
- With the updated particle locations, the blob solution is constructed

```matlab
for i = 1:Nr
    f(i) = sum(w.*psi_1d(xr(i)-x,epsilon));
end
```
or 
```matlab
parfor i = 1:Nr
    for j = 1:Nr
        f(i,j) = sum(W.*psi_2d(vrx(i,j)-Vx,vry(i,j)-Vy,epsilon));
    end
end
```
Using the particle locations and particle weights, the blob solution at the current time is computed and overwritten into `f`

For the 1D examples 
```matlab
for i = 1:Nr
    f(i) = sum(w.*psi_1d(xr(i)-x,epsilon));
end
```
and for the 2D examples
```matlab
parfor i = 1:Nr
    for j = 1:Nr
        f(i,j) = sum(W.*psi_2d(vrx(i,j)-Vx,vry(i,j)-Vy,epsilon));
    end
end
```

#### The `exact` and `exact_2d` functions

- The arguments for the `exact` function are:
    - `x` is an array of locations where the exact solution is sampled at
    - `t` is the current time
    - `m` is an argument if you are working with the porous medium example
- The output of `exact` is
    - `f` an array of the exact solution sampled at `t` and `x`.

- The arguments for the `exact_2d` function are:
  - `vx` is an array of $x$ coordinate locations where the exact solution is sampled at
  - `vy` is an array of $y$ coordinate locations where the exact solution is sampled at
  - `t` is the current time
- The output of `exact_2d` is
  - `f` is an array of the exact solution sampled at `t` and `vx` and `vy`.

- The energy dissipation, fixed point iteration information, and any important conserved quantities are computed and saved to error_list.  The computation of the errors rely on the `exact` or `exact_2d` functions to compute the exact function sampled on the reference mesh.  In one dimension this is done with the following code:

```matlab
f_exact = exact(time,xr);
error_list(nt,1) = time;
f_error = f-f_exact;

%relative error
Linf_error = max(abs(f_error))/max(abs(f_exact));
L1_error = sum(abs(f_error))/sum(abs(f_exact));
L2_error = sqrt(sum(f_error.^2)/sum(f_exact.^2));
error_list(nt,2) = Linf_error;
error_list(nt,3) = L1_error;
error_list(nt,4) = L2_error;

rho = sum(w);
inside = zeros(1,n);
for i = 1:n
    inside(i) = sum(w.*psi_1d(xr(i)-x,epsilon));
end
eta = dx*sum(inside.*log(inside));
error_list(nt,5) = rho;
error_list(nt,6) = eta;
```

In the 2D BKW case, the errors are computed since there is an exact solution to compare with.  In the Coulomb case the error is not computed.
```matlab
f_error = f-f_exact;
% relative error
Linf_error = max(max(abs(f_error)))/max(max(abs(f_exact)));
L1_error = sum(sum(abs(f_error)))/sum(sum(abs(f_exact)));
L2_error = sqrt(sum(sum(f_error.^2))/sum(sum(f_exact.^2)));   
error_list(nt,2) = Linf_error;   
error_list(nt,3) = L1_error;   
error_list(nt,4) = L2_error;

% moments
rho = sum(W);
m1 = sum(W.*Vx);
m2 = sum(W.*Vy);
E = sum(W.*(Vx.^2+Vy.^2));   
F = zeros(Nr^2,1);
% parfor for i
parfor i = 1:Nr^2
    F(i) = sum(W.*psi_2d(Vrx(i)-Vx,Vry(i)-Vy,epsilon));
end
%compute the fisher information 
term1_x = zeros(Np,1);
term1_y = zeros(Np,1);
for i = 1:Np %compute fisher information and dissipation terms
    % at time t = t0 + nt*dt
    [A,B] = gpsi_2d(Vx(i)-Vrx,Vy(i)-Vry,epsilon);
    term1_x(i) = dv^2*sum(A.*log(F));
    term1_y(i) = dv^2*sum(B.*log(F));

end
Fish = sum(W.*(term1_x.^2 + term1_y.^2)); %fisher information 
eta = dv^2*sum(F.*log(F));
eta1 = dv^2*sum(F.*(log(F)+log(2*pi)+(Vrx.^2+Vry.^2)/2));
error_list(nt,5) = rho;   
error_list(nt,6) = m1;   
error_list(nt,7) = m2;
error_list(nt,8) = E;  
error_list(nt,9) = eta;  
error_list(nt,10) = eta1;
error_list(nt,11) = Fish;
error_list(nt,12) = dissipation;
```
## Reproducing Results in the Paper

To reproduce the results given in section 4 of the paper, open and run the 
`particle_method` or `particle_method_2d_parallel` in the folder corresponding 
to the example for which you would like to reproduce results (see contents listed above).  For each entry of `n_list` 
defined at the beginning of `particle_method` or `particle_method_2d_parallel`, a `.mat` file will be saved into the `Data` 
directory located in that folder.  Next, open and run the `plots` file to create the plots given in the paper.


If you would like to alter the parameters listed outlined in the project description section above, run the `particle_method` 
`particle_method_2d_parallel` files which will save `.mat` files to the data `Data` repository.  Depending on what parameters are altered, the title of the `.mat` file that is saved in `Data` may not match the title that is loaded in `plots` and thus you may need 
to alter the argument of `load()` to match this saved file.



