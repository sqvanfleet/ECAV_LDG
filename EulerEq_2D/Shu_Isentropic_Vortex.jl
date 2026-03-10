using StartUpDG
using LinearAlgebra 
using StaticArrays
using OrdinaryDiffEqTsit5, RecursiveArrayTools
using JLD2
using Trixi
using Trixi.ForwardDiff
using Plots
using JLD2
using Printf

initial_condition_type = :shu_isentropic_vortex

include("initial_conditions.jl")
include("rhs.jl")

psi_1(u, normal, ::CompressibleEulerEquations2D) = u[2] * normal
psi_2(u, normal, ::CompressibleEulerEquations2D) = u[3] * normal

dudv(v, equations) = ForwardDiff.jacobian(v -> entropy2cons(v, equations), v)

gamma = 1.4
equations = CompressibleEulerEquations2D(gamma)


#shu_isentropic_euler_vortex Problem
if initial_condition_type == :shu_isentropic_vortex

        N_list = [1,2,3,4]
        K_list = [16,32,64,128]

        DG_type = :modal
        AV_discretization = :LDG

        for N in N_list
                for K1D in K_list
                        #Mesh Parameters
                        # N = 1 #Order of Polynomial
                        # K1D = 8 #Number of elements in each direction
                        display("N = $(N), K = $(K1D)")
                        if DG_type == :nodal
                                local rd = RefElemData(Tri(), SBP(), N)
                        elseif DG_type == :modal
                                local rd = RefElemData(Tri(), N)
                        end

                        local VXY, EToV = uniform_mesh(rd.element_type, K1D)
                        local L = 10.0 #Length of the domain
                        local VXY = VXY.*L #Scale the mesh
                        local md = MeshData((VXY,EToV), rd; is_periodic = true)
                        local x,y = md.x, md.y
                        local u = rd.Pq * shu_isentropic_vortex.(md.xq,md.yq,0.0,L,gamma,equations)
                        local Final_time = 15.0
                        #local Final_time = sqrt(2)*5.0/sqrt(2/gamma)
                        #local Final_time = 2.0*sqrt(gamma)*L - 5.0
                        local tspan = (0.0, Final_time)
                        local du = similar(VectorOfArray(u))
                        local params = (; rd, md, equations, interface_flux = flux_hllc, AV_discretization, gamma)
                        local ode = ODEProblem(rhs!, VectorOfArray(u), tspan, params)
                        local sol = solve(ode, Tsit5(); abstol = 1e-10, reltol = 1e-8, 
                                        saveat=LinRange(tspan..., 50), 
                                                callback=AliveCallback(alive_interval=10))
                        local u = parent(sol.u[end])


                        # higher accuracy quadrature rule
                        local rq, sq, wq = quad_nodes(rd.element_type, rd.N + 5)
                        local Vq = vandermonde(rd.element_type, rd.N, rq, sq) / rd.VDM
                        local xq, yq = Vq * md.x, Vq * md.y

                        # local xp, yp, up = rd.Vp * x, rd.Vp * y, rd.Vp * getindex.(u,1)

                        # display(scatter(vec(md.x), vec(md.y), zcolor=vec(getindex.(u,1)), msw=0, ms = 1, legend=false, ratio=1, cam=(0,90),
                        #         title="DG_Density", xlabel="x", ylabel="y", c=:red))

                        local alpha = pi/4
                        local M_inf = sqrt(2.0/gamma)
                        local sigma = 1.0
                        local R = 1.0
                        local Beta = M_inf*5.0*sqrt(2.0)*exp(0.5)/(4*pi)
                        local ux = M_inf * cos(alpha)
                        local vx = M_inf * sin(alpha)

                        local u_exact = shu_isentropic_vortex.(xq .- ux*sol.t[end], yq .- vx*sol.t[end], 0.0, L, gamma, equations)
                        # local u_exact = shu_isentropic_vortex.(xq, yq, sol.t[end], L, gamma, equations)
                        local error = Vq * u - u_exact
                        # up_exact = rd.Vp * getindex.(u_exact,1)

                        # display(scatter(vec(xp), vec(yp), zcolor=vec(up_exact), msw=0, ms = 1, legend=false, ratio=1, cam=(0,90),
                        # title="exact Density", xlabel="x", ylabel="y", c=:red))

                        # L2_norm = sqrt(sum(md.wJq .* norm.(rd.Vq * u - u_exact)).^2)
                        local wJq = wq.*Vq*md.J
                        local L2_norm = sqrt(sum(wJq .* norm.(error).^2))

                        local L2_rel_norm = L2_norm / sqrt(sum(wJq .* norm.(u_exact).^2))

                        display("N = $(N), L2 norm = $(L2_norm)")

                        local base_path = @__DIR__

                        local folder = joinpath(base_path, "Vortex_Data_10")
                        mkpath(folder)
                        local filename = "N$(N)_K$(K1D)_T$(sol.t[end])_$(DG_type)_$(AV_discretization)_$(sol.prob.p.interface_flux).jld2"
                        
                        @save joinpath(folder, filename) md rd sol u_exact L2_norm L2_rel_norm
                        
                end
        end
end



#Enter yes if you want to generate a convergence table.
ConvergenceTable = :yes

if ConvergenceTable == :yes
    #Convergence Table
    # Define the values of K
    K_vals = [16,32,64,128]

    # Helper function to format in LaTeX scientific notation
    function latex_sci_notation(x::Float64; sig_digits=4)
        # Create format string, e.g. "%.4e"
        fmt = "%" * "." * string(sig_digits - 1) * "e"
        str = Printf.format(Printf.Format(fmt), x)
        base, exp = split(str, 'e')
        exp_val = parse(Int, exp)
        return "\$$(base) \\times 10^{ $(exp_val) }\$"
    end

    function generate_latex_table_multi(K_vals, data_series)
        num_orders = length(data_series)
        header_cols = ["\$K\$"]
        for data in data_series
            push!(header_cols, "\$L^2\$ (" * data.label * ")")
            push!(header_cols, "Order")
        end

        latex_table = "\\begin{tabular}{|" * join(fill("c", length(header_cols)), "|") * "|}\n"
        latex_table *= "\\hline\n"
        latex_table *= join(header_cols, " & ") * " \\\\\n"
        latex_table *= "\\hline\n"

        num_rows = length(K_vals)

        for i in 1:num_rows
            row = ["\$$(K_vals[i])\$"]
            for data in data_series
                L2 = data.L2_error[i]
                err_str = latex_sci_notation(L2)
                if i == 1
                    order_str = "---"
                else
                    order = data.Approx_Order[i - 1]
                    order_str = "\$" * @sprintf("%.4f", order) * "\$"
                end
                push!(row, err_str)
                push!(row, order_str)
            end
            latex_table *= join(row, " & ") * " \\\\\n"
        end

        latex_table *= "\\hline\n\\end{tabular}"
        return latex_table
    end

    function generate_latex_table_stacked(K_vals, data_series)
        num_K = length(K_vals)
        num_series = length(data_series)
        half = num_series ÷ 2  # Assume even number of data_series

        # Build header
        header_cols_top = ["\$K\$"]
        for data in data_series[1:half]
            push!(header_cols_top, "\$L^2\$ (" * data.label * ")")
            push!(header_cols_top, "Order")
        end

        header_cols_bottom = ["\$K\$"]
        for data in data_series[half+1:end]
            push!(header_cols_bottom, "\$L^2\$ (" * data.label * ")")
            push!(header_cols_bottom, "Order")
        end

        # Begin LaTeX table
        num_cols = 1 + (half * 2)
        latex_table = "\\begin{tabular}{|" * join(fill("c", num_cols), "|") * "|}\n"
        latex_table *= "\\hline\n"
        latex_table *= join(header_cols_top, " & ") * " \\\\\n"
        latex_table *= "\\hline\n"
 
        # Top half rows
        for i in 1:num_K
            row = ["\$$(K_vals[i])\$"]
            for data in data_series[1:half]
                L2 = data.L2_error[i]
                err_str = latex_sci_notation(L2)
                if i == 1
                    order_str = "---"
                else
                    order = data.Approx_Order[i - 1]
                    order_str = "\$" * @sprintf("%.4f", order) * "\$"
                end
                push!(row, err_str)
                push!(row, order_str)
            end
            latex_table *= join(row, " & ") * " \\\\\n"
        end

        # Bottom half header row
        latex_table *= "\\hline\n"
        latex_table *= join(header_cols_bottom, " & ") * " \\\\\n"
        latex_table *= "\\hline\n"

        # Bottom half rows
        for i in 1:num_K
            row = ["\$$(K_vals[i])\$"]
            for data in data_series[half+1:end]
                L2 = data.L2_error[i]
                err_str = latex_sci_notation(L2)
                if i == 1
                    order_str = "---"
                else
                    order = data.Approx_Order[i - 1]
                    order_str = "\$" * @sprintf("%.4f", order) * "\$"
                end
                push!(row, err_str)
                push!(row, order_str)
            end
            latex_table *= join(row, " & ") * " \\\\\n"
        end

        latex_table *= "\\hline\n\\end{tabular}"
        return latex_table
    end

    #Order 1

    data_path = joinpath(@__DIR__, "Vortex_Data_10")

    L2_error_K16_N1 = load(joinpath(data_path, "N1_K16_T15.0_$(DG_type)_LDG_flux_hllc.jld2"), "L2_rel_norm");
    L2_error_K32_N1 = load(joinpath(data_path, "N1_K32_T15.0_$(DG_type)_LDG_flux_hllc.jld2"), "L2_rel_norm");
    L2_error_K64_N1 = load(joinpath(data_path, "N1_K64_T15.0_$(DG_type)_LDG_flux_hllc.jld2"), "L2_rel_norm");
    L2_error_K128_N1 = load(joinpath(data_path, "N1_K128_T15.0_$(DG_type)_LDG_flux_hllc.jld2"), "L2_rel_norm");
    
    error_ratio_N1 = [
        L2_error_K16_N1 / L2_error_K32_N1,
        L2_error_K32_N1 / L2_error_K64_N1,
        L2_error_K64_N1 / L2_error_K128_N1]

    Approx_Order_N1 = log.(error_ratio_N1) ./ log(2)

    L2_error_N1 = [
        L2_error_K16_N1,
        L2_error_K32_N1,
        L2_error_K64_N1,
        L2_error_K128_N1]

    L2_error_K16_N2 = load(joinpath(data_path, "N2_K16_T15.0_$(DG_type)_LDG_flux_hllc.jld2"), "L2_rel_norm");
    L2_error_K32_N2 = load(joinpath(data_path, "N2_K32_T15.0_$(DG_type)_LDG_flux_hllc.jld2"), "L2_rel_norm");
    L2_error_K64_N2 = load(joinpath(data_path, "N2_K64_T15.0_$(DG_type)_LDG_flux_hllc.jld2"), "L2_rel_norm");
    L2_error_K128_N2 = load(joinpath(data_path, "N2_K128_T15.0_$(DG_type)_LDG_flux_hllc.jld2"), "L2_rel_norm");
   
    error_ratio_N2 = [
        L2_error_K16_N2 / L2_error_K32_N2,
        L2_error_K32_N2 / L2_error_K64_N2,
        L2_error_K64_N2 / L2_error_K128_N2]

    Approx_Order_N2 = log.(error_ratio_N2) ./ log(2)

    L2_error_N2 = [
        L2_error_K16_N2,
        L2_error_K32_N2,
        L2_error_K64_N2,
        L2_error_K128_N2]


    #Order =3

    L2_error_K16_N3 = load(joinpath(data_path, "N3_K16_T15.0_$(DG_type)_LDG_flux_hllc.jld2"), "L2_rel_norm");
    L2_error_K32_N3 = load(joinpath(data_path, "N3_K32_T15.0_$(DG_type)_LDG_flux_hllc.jld2"), "L2_rel_norm");
    L2_error_K64_N3 = load(joinpath(data_path, "N3_K64_T15.0_$(DG_type)_LDG_flux_hllc.jld2"), "L2_rel_norm");
    L2_error_K128_N3 = load(joinpath(data_path, "N3_K128_T15.0_$(DG_type)_LDG_flux_hllc.jld2"), "L2_rel_norm");     


    error_ratio_N3 = [
        L2_error_K16_N3 / L2_error_K32_N3,
        L2_error_K32_N3 / L2_error_K64_N3,
        L2_error_K64_N3 / L2_error_K128_N3]

    Approx_Order_N3 = log.(error_ratio_N3) ./ log(2)

    L2_error_N3 = [
        L2_error_K16_N3,
        L2_error_K32_N3,
        L2_error_K64_N3,
        L2_error_K128_N3]

    #Order 4
    L2_error_K16_N4 = load(joinpath(data_path, "N4_K16_T15.0_$(DG_type)_LDG_flux_hllc.jld2"), "L2_rel_norm");
    L2_error_K32_N4 = load(joinpath(data_path, "N4_K32_T15.0_$(DG_type)_LDG_flux_hllc.jld2"), "L2_rel_norm");
    L2_error_K64_N4 = load(joinpath(data_path, "N4_K64_T15.0_$(DG_type)_LDG_flux_hllc.jld2"), "L2_rel_norm");
    L2_error_K128_N4 = load(joinpath(data_path, "N4_K128_T15.0_$(DG_type)_LDG_flux_hllc.jld2"), "L2_rel_norm");

    error_ratio_N4 = [
        L2_error_K16_N4 / L2_error_K32_N4,
        L2_error_K32_N4 / L2_error_K64_N4,
        L2_error_K64_N4 / L2_error_K128_N4]

    Approx_Order_N4 = log.(error_ratio_N4) ./ log(2)

    L2_error_N4 = [
        L2_error_K16_N4,
        L2_error_K32_N4,
        L2_error_K64_N4,
        L2_error_K128_N4]

    data_N1 = (L2_error=L2_error_N1, Approx_Order=Approx_Order_N1, label="N=1")
    data_N2 = (L2_error=L2_error_N2, Approx_Order=Approx_Order_N2, label="N=2")
    data_N3 = (L2_error=L2_error_N3, Approx_Order=Approx_Order_N3, label="N=3")
    data_N4 = (L2_error=L2_error_N4, Approx_Order=Approx_Order_N4, label="N=4")

    data_series = [data_N1 data_N2 data_N3 data_N4]

    #latex_table = generate_latex_table_multi(K_vals, data_series)
    latex_table = generate_latex_table_stacked(K_vals,data_series)
    println(latex_table) 
end




