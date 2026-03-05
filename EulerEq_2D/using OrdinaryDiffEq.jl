using OrdinaryDiffEq

f(u,p,t) = -u
u0 = [1.0, 2.0, 3.0]
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)

sol = solve(prob, RDPK3Sp510(); abstol=1e-8, reltol=1e-6,
            saveat=LinRange(tspan..., 10),
            callback=AliveCallback(alive_interval=10))

println(sol[end])