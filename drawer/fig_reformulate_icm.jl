using CSV, DataFrames
using CairoMakie, LaTeXStrings

begin
    data_file = "data/ICM_Ewald3D_2-1.csv"
    data_file_force = "data/reformulate_ewald_force.csv"

    df = CSV.read(data_file, DataFrame)
    df_force = CSV.read(data_file_force, DataFrame)

    N_imgs = unique(df.N_img)
    N_slabs = unique(df.N_slab)

    f = Figure(backgroundcolor = RGBf(1.0, 1.0, 1.0), size = (500, 700), fontsize = 20)
    ga = f[1, 1]
    gb = f[2, 1]
    axr = Axis(ga, xlabel = L"M", ylabel = L"$\mathcal{E}_r$ (Energy)", title = L"\gamma_u = \gamma_d = 0.95", yscale = log10)
    axl = Axis(gb, xlabel = L"M", ylabel = L"$\mathcal{E}_r$ (Force)", yscale = log10)
    ylims!(axr, 1e-11, 1e2)
    ylims!(axl, 1e-11, 1e2)

    colors = [:green, :blue, :red]

    Ns = [4, 7, 10]
    H = 5.0
    Δ = 0.95
    γ = (Δ, Δ)
    for (i, N_slab) in enumerate(N_slabs)
        Lz = H * (2 * Ns[i] + 1)
        dfγ = df[(df.N_slab .== N_slab) .& (df.γ1 .== γ[1]) .& (df.γ2 .== γ[2]), :]
        scatter!(axr, dfγ.N_img, dfγ.error_r, label = L"L_z/H = %$(2 * Ns[i] + 1)", color = colors[i])

        p = [0.01, 1.0]

        fe = x -> 
            p[1] * exp(-(Lz - H) * 2π / 15) + 
            # p[2] * exp(-(Lz - H)^2) + 
            p[1] * sum([2 *  Δ^(j) * exp(-(Lz - (j + 1) * H) * 2π / 15) for j in 1:Int(x)]) + 
            # p[2] * sum([2 * Δ^(j) * exp(-(Lz - (j + 1) * H)^2) for j in 1:Int(x)]) + 
            p[2] * Δ^(x + 1) * exp(-2π * (x + 1) * 5 / 15)
            
        lines!(axr, [2:16...], fe.([2:16...]), color = colors[i], linestyle = :dash)

        df_forceγ = df_force[(df_force.N_slab .== N_slab) .& (df_force.γ1 .== γ[1]) .& (df_force.γ2 .== γ[2]), :]
        scatter!(axl, df_forceγ.N_img, df_forceγ.error_r, color = colors[i])
        
        p = [0.2, 3.0]
        lines!(axl, [2:16...], fe.([2:16...]), color = colors[i], linestyle = :dash)
    end

    axislegend(axr, position = :lt, labelsize = 18, merge = true, unique = true)

    text!(axr, (2, 10^(-9)), text = "(a)", fontsize = 30, align = (:left, :baseline))
    text!(axl, (2, 10^(-9)), text = "(b)", fontsize = 30, align = (:left, :baseline))

    save("figs/reformulate_icm_ewald.png", f)
    save("figs/reformulate_icm_ewald.pdf", f)

    f
end