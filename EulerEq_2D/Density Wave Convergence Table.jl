using JLD2
using StartUpDG
using LinearAlgebra 
using StaticArrays
using Plots
using Printf
using Trixi

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
    K_vals = [4, 8, 16, 32, 64]

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

    #Order 1
    L2_norm_K4_N1 = load("Density_wave_Data/N1_K4_$(DG_type)_$(AV_discretization).jld2", "L2_rel_norm");
    L2_norm_K8_N1 = load("Density_wave_Data/N1_K8_$(DG_type)_$(AV_discretization).jld2", "L2_rel_norm");
    L2_norm_K16_N1 = load("Density_wave_Data/N1_K16_$(DG_type)_$(AV_discretization).jld2", "L2_rel_norm");
    L2_norm_K32_N1 = load("Density_wave_Data/N1_K32_$(DG_type)_$(AV_discretization).jld2", "L2_rel_norm");
    L2_norm_K64_N1 = load("Density_wave_Data/N1_K64_$(DG_type)_$(AV_discretization).jld2", "L2_rel_norm");

    error_ratio_N1 = [
        L2_norm_K4_N1 / L2_norm_K8_N1,
        L2_norm_K8_N1 / L2_norm_K16_N1,
        L2_norm_K16_N1 / L2_norm_K32_N1,
        L2_norm_K32_N1 / L2_norm_K64_N1]

    Approx_Order_N1 = log.(error_ratio_N1) ./ log(2)

    L2_error_N1 = [
        L2_norm_K4_N1,
        L2_norm_K8_N1,
        L2_norm_K16_N1,
        L2_norm_K32_N1,
        L2_norm_K64_N1]

    # Order 2
    L2_norm_K4_N2 = load("Density_wave_Data/N2_K4_$(DG_type)_$(AV_discretization).jld2", "L2_norm");
    L2_norm_K8_N2 = load("Density_wave_Data/N2_K8_$(DG_type)_$(AV_discretization).jld2", "L2_norm");
    L2_norm_K16_N2 = load("Density_wave_Data/N2_K16_$(DG_type)_$(AV_discretization).jld2", "L2_norm");
    L2_norm_K32_N2 = load("Density_wave_Data/N2_K32_$(DG_type)_$(AV_discretization).jld2", "L2_norm");
    L2_norm_K64_N2 = load("Density_wave_Data/N2_K64_$(DG_type)_$(AV_discretization).jld2", "L2_norm");


    error_ratio_N2 = [
        L2_norm_K4_N2 / L2_norm_K8_N2,
        L2_norm_K8_N2 / L2_norm_K16_N2,
        L2_norm_K16_N2 / L2_norm_K32_N2,
        L2_norm_K32_N2 / L2_norm_K64_N2]

    Approx_Order_N2 = log.(error_ratio_N2) ./ log(2)

    L2_error_N2 = [
        L2_norm_K4_N2,
        L2_norm_K8_N2,
        L2_norm_K16_N2,
        L2_norm_K32_N2,
        L2_norm_K64_N2]
    # #Order =3
    # L2_norm_K2_N3 = load("Density_wave_Data/N3_K2_$(DG_type)_$(AV_discretization).jld2", "L2_norm");
    # L2_norm_K4_N3 = load("Density_wave_Data/N3_K2_$(DG_type)_$(AV_discretization).jld2", "L2_norm");
    # L2_norm_K8_N3 = load("Density_wave_Data/N3_K8_$(DG_type)_$(AV_discretization).jld2", "L2_norm");
    # L2_norm_K16_N3 = load("Density_wave_Data/N3_K16_$(DG_type)_$(AV_discretization).jld2", "L2_norm");
    # L2_norm_K32_N3 = load("Density_wave_Data/N3_K32_$(DG_type)_$(AV_discretization).jld2", "L2_norm");


    # error_ratio_N3 = [
    #     L2_norm_K2_N3 / L2_norm_K4_N3,
    #     L2_norm_K4_N3 / L2_norm_K8_N3,
    #     L2_norm_K8_N3 / L2_norm_K16_N3,
    #     L2_norm_K16_N3 / L2_norm_K32_N3]

    # Approx_Order_N3 = log.(error_ratio_N3) ./ log(2)

    # L2_error_N3 = [
    #     L2_norm_K2_N3,
    #     L2_norm_K4_N3,
    #     L2_norm_K8_N3,
    #     L2_norm_K16_N3,
    #     L2_norm_K32_N3]

    # #Order 4
    # L2_norm_K2_N4 = load("Density_wave_Data/N4_K4_$(DG_type)_$(AV_discretization).jld2", "L2_norm");
    # L2_norm_K4_N4 = load("Density_wave_Data/N4_K4_$(DG_type)_$(AV_discretization).jld2", "L2_norm");
    # L2_norm_K8_N4 = load("Density_wave_Data/N4_K8_$(DG_type)_$(AV_discretization).jld2", "L2_norm");
    # L2_norm_K16_N4 = load("Density_wave_Data/N4_K16_$(DG_type)_$(AV_discretization).jld2", "L2_norm");
    # L2_norm_K32_N4 = load("Density_wave_Data/N4_K32_$(DG_type)_$(AV_discretization).jld2", "L2_norm");


    # error_ratio_N4 = [
    #     L2_norm_K2_N4 / L2_norm_K4_N4,
    #     L2_norm_K4_N4 / L2_norm_K8_N4,
    #     L2_norm_K8_N4 / L2_norm_K16_N4,
    #     L2_norm_K16_N4 / L2_norm_K32_N4]

    # Approx_Order_N4 = log.(error_ratio_N4) ./ log(2)

    # L2_error_N4 = [
    #     L2_norm_K2_N4,
    #     L2_norm_K4_N4,
    #     L2_norm_K8_N4,
    #     L2_norm_K16_N4,
    #     L2_norm_K32_N4]

    data_N1 = (L2_error=L2_error_N1, Approx_Order=Approx_Order_N1, label="N=1")
    data_N2 = (L2_error=L2_error_N2, Approx_Order=Approx_Order_N2, label="N=2")
    # data_N3 = (L2_error=L2_error_N3, Approx_Order=Approx_Order_N3, label="N=3")
    # data_N4 = (L2_error=L2_error_N4, Approx_Order=Approx_Order_N4, label="N=4")

    data_series = [data_N1, data_N2] 

    # data_series = [data_N1, data_N2, data_N3, data_N4]

    latex_table = generate_latex_table_multi(K_vals, data_series)
    println(latex_table)
end

# # Plot the error, entropy, epsilon evolution
# N = 2
# K = 32

# #BR1 versus LDG switch
# #AV_discretization = :BR1
# AV_discretization = :LDG

# #Choose Modal or Nodal
# DG_type = :nodal
# #DG_type = :modal

# if DG_type == :nodal && AV_discretization == :LDG
#     #Load the data for nodal LDG
#     md, rd, sol = load("Density_wave_Data/N$(N)_K$(K1D)_nodal_LDG.jld2", "md", "rd", "sol")
# elseif DG_type == :nodal && AV_discretization == :BR1
#     #Load the data for nodal BR1
#     md, rd, sol = load("Density_wave_Data/N$(N)_K$(K1D)_nodal_BR1.jld2", "md", "rd", "sol")
# elseif DG_type == :modal && AV_discretization == :LDG
#     #Load the data for modal LDG
#     md, rd, sol = load("Density_wave_Data/N$(N)_K$(K1D)_modal_LDG.jld2", "md", "rd", "sol")
# elseif DG_type == :modal && AV_discretization == :BR1
#     #Load the data for modal BR1
#     md, rd, sol = load("Density_wave_Data/N$(N)_K$(K1D)_modal_BR1.jld2", "md", "rd", "sol")
# else            
#     error("Invalid combination of DG_type and AV_discretization")
# end

# function density_wave(x,y,equations::CompressibleEulerEquations2D)
#     rho = 1.0 + 0.5 * sin(2 * pi * (x + y))
#     u = 0.1
#     v = 0.2
#     p = 10.0
#     return prim2cons(SVector(rho, u, v, p), equations)
# end

# x, y = md.x, md.y

# error_evolution = zeros(length(sol.t))






