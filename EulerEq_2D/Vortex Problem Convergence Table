using JLD2
# using StartUpDG
# using LinearAlgebra 
# using StaticArrays
# using Plots
using Printf
using Trixi
using Trixi.ForwardDiff
using OrdinaryDiffEq, RecursiveArrayTools
using LaTeXStrings

gamma = 1.4
equations = CompressibleEulerEquations2D(gamma)

AV_discretization = :LDG
#AV_discretization = :BR1



#DG_type = :nodal
DG_type = :modal

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




function periodic_wrap(x, a, T)
    return mod(x - a, T) + a
end

function shu_isentropic_euler_vortex(x,y,t,L,gamma,equations::CompressibleEulerEquations2D)
    alpha = pi/4
    M_inf = sqrt(2.0/gamma)
    sigma = 1.0
    R = 1.0
    Beta = M_inf*5.0*sqrt(2.0)*exp(0.5)/(4*pi)
    ux = M_inf * cos(alpha)
    vx = M_inf * sin(alpha)
    #t = periodic_wrap(t, 0.0, 2.0*L*sqrt(gamma))
    x = periodic_wrap(x, -L, 2.0*L)
    y = periodic_wrap(y, -L, 2.0*L)
    r = sqrt(((x/R)-ux*t)^2 + ((y/R)-vx*t)^2)
    Omega(r) = Beta * exp(-r^2/(2.0*sigma^2)) 
    rho = (1.0 - 0.5 * (gamma - 1.0) * Omega(r)^2)^(1.0/(gamma - 1.0))
    u = M_inf * cos(alpha) - (y/R - vx*t)*Omega(r)
    v = M_inf * sin(alpha) + (x/R - ux*t)*Omega(r)
    p = (1/gamma)*(1.0 - 0.5 * (gamma - 1.0) * Omega(r)^2)^(gamma/(gamma - 1.0))
    return prim2cons(SVector(rho, u, v, p), equations)
end
