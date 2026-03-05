using StartUpDG
using LinearAlgebra 
using StaticArrays
using OrdinaryDiffEq, RecursiveArrayTools
using Trixi
using Trixi.ForwardDiff
 
function periodic_wrap(x, a, T)
    return mod(x - a, T) + a
end

function Shock_vortex_interaction(x,y,equations::CompressibleEulerEquations2D,gamma)
        rho_L = 1.0
        u_L = sqrt(gamma)
        v_L = 0.0
        p_L = 1.0
        Ms = 1.1

        RH1 = (2.0 + (gamma - 1.0) * Ms^2)/ ((gamma +1.0)*Ms^2)
        RH2 = 1.0 + ((2.0*gamma)/(gamma + 1.0)) * (Ms^2 - 1.0)

        rho_R = (1/RH1)*rho_L
        u_R = RH1*u_L
        v_R = 0.0
        p_R = RH2*p_L

        x_c = 0.25
        y_c = 0.5
        r_c = 0.05
        alpha = 0.204
        epsilon = 0.3

        r = sqrt((x - x_c)^2 + (y - y_c)^2)
        theta = atan((y - y_c), (x - x_c))
        tao = r/r_c

        v_theta = epsilon*tao*exp(alpha*(1.0 - tao^2))
        delta_u = v_theta * sin(theta)
        delta_v = -v_theta * cos(theta)
        delta_T = -((gamma - 1.0)*epsilon^2*exp(2.0*alpha*(1.0 - tao^2)))/(4.0*alpha*gamma)
        T_L = p_L/rho_L
        T_vor = delta_T + T_L

        if x < 0.5
                # rho = rho_L
                # u = u_L
                # v = v_L
                # p = p_L
                rho = rho_L*(T_vor/T_L)^(1.0/(gamma - 1.0))
                u = u_L + delta_u
                v = v_L + delta_v
                p = p_L*(T_vor/T_L)^(1.0/(gamma - 1.0))
        else x > 0.5
                # rho = rho_R
                # u = u_R
                # v = v_R
                # p = p_R
                rho = rho_R*(T_vor/T_L)^(1.0/(gamma - 1.0))
                u = u_R + delta_u
                v = v_R + delta_v
                p = p_R*(T_vor/T_L)^(1.0/(gamma - 1.0))
        end
        return prim2cons(SVector(rho, u, v, p),equations)
end

function shu_isentropic_vortex(x,y,t,L,gamma,equations::CompressibleEulerEquations2D)
    alpha = pi/4
    M_inf = sqrt(2.0/gamma)
    sigma = 1.0
    R = 1.0
    Beta = M_inf*5.0*sqrt(2.0)*exp(0.5)/(4*pi)
    ux = M_inf * cos(alpha)
    vx = M_inf * sin(alpha)
    #t = periodic_wrap(t, 0.0, 2.0*L*sqrt(gamma))
    x = periodic_wrap(x-ux*t, -L, 2.0*L)
    y = periodic_wrap(y-vx*t, -L, 2.0*L)
    r = sqrt(((x/R))^2 + ((y/R))^2)
    Omega(r) = Beta * exp(-r^2/(2.0*sigma^2)) 
    rho = (1.0 - 0.5 * (gamma - 1.0) * Omega(r)^2)^(1.0/(gamma - 1.0))
    u = M_inf * cos(alpha) - (y/R)*Omega(r)
    v = M_inf * sin(alpha) + (x/R)*Omega(r)
    p = (1/gamma)*(1.0 - 0.5 * (gamma - 1.0) * Omega(r)^2)^(gamma/(gamma - 1.0))
    return prim2cons(SVector(rho, u, v, p), equations)
end

function density_wave(x,y,t,equations::CompressibleEulerEquations2D)
    
    u = 0.1
    v = 0.2
    rho = 1.0 + 0.5 * sin(2 * pi * ((x-u*t) + (y-v*t)))
    p = 10.0
    return prim2cons(SVector(rho, u, v, p), equations)
end

function Riemann_problem_config_1(x,y,equations::CompressibleEulerEquations2D)
        if x>0 && y>0
                rho, u, v, p = 0.5313, 0.0, 0.0, 0.4
        elseif x<0 && y>0
                rho, u, v, p = 1.0, 0.7276, 0.0, 1.0
        elseif x<0 && y<0
                rho, u, v, p = 0.8, 0.0, 0.0, 1.0
        elseif x>0 && y<0
                rho, u, v, p = 1.0, 0.0, 0.7276, 1.0
        end
        return prim2cons(SVector(rho, u, v, p), equations)
end