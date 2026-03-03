using Trixi

function shu_osher(x, equations::CompressibleEulerEquations1D)
    if x < -4
        return prim2cons(SVector(3.857143, 2.629369, 10.33333), equations)
    else
        rho = 1 + 0.2 * sin(5 * x)
        u = 0.0
        p = 1.0
        return prim2cons(SVector(rho, u, p), equations)
    end
end

function Euler_problem_1(x, equations::CompressibleEulerEquations1D, k, pmin, rhomin)
    rho = rhomin + exp(2.0*sin(1.0 + k*pi*x))
    u = 0.1*cos(1.0 + k*pi*x)
    p = pmin + 0.5*(1.0 - cos(k*pi*x - 0.25))
    return prim2cons(SVector(rho, u, p), equations)
end

function stationary_contact_wave(x, equations::CompressibleEulerEquations1D; 
                                        amplitude = 0.5)
    rho = 1.0 + amplitude * (sin(2*pi*x) + (abs(x) < 0.3))
    # if abs(x) < 0.3
    #     rho = 1.0 + amplitude
    # else
    #     rho = 1.0
    # end
    u = 0.0
    p = 1.0

    return prim2cons(SVector(rho, u, p), equations)
end

function initial_condition_density_wave(x, equations::CompressibleEulerEquations1D; 
                                        amplitude = 0.5, a = 10.0)
    rho = 1.0 + 0.5 * exp(-a * sin(pi * x)^2)
    # rho = 1.0 + amplitude * sin(20*pi*x)
    # rho = 1.0 + amplitude * (sin(20*pi*x) + (abs(x) < 0.3))

    u = 0.1
    p = 10.0
    return prim2cons(SVector(rho, u, p), equations)
end

function density_wave_convergence(x, equations::CompressibleEulerEquations1D)
    rho = 1.0 + 0.5 * sin(2.0 * pi * x)
    u = 0.1
    p = 10.0
    return prim2cons(SVector(rho, u, p), equations)
end

function initial_condition_sod_shock_tube(x, equations::CompressibleEulerEquations1D)
    if x < 0.3
        return prim2cons(SVector(1.0, 0.75, 1.0), equations)
    else
        return prim2cons(SVector(0.125, 0.0, 0.1), equations)
    end
end