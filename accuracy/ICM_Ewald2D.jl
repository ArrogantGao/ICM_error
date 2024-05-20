using ExTinyMD, EwaldSummations, CairoMakie, LaTeXStrings

n_atoms = 32
Lxy = 10.0
Lz = 1.0
L = (Lxy, Lxy, Lz)
boundary = ExTinyMD.Q2dBoundary(Lxy, Lxy, Lz)

atoms = Vector{Atom{Float64}}()
for i in 1:n_atoms÷2
    push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
end

for i in n_atoms÷2 + 1 : n_atoms
    push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
end

info = SimulationInfo(n_atoms, atoms, (0.0, Lxy, 0.0, Lxy, 0.0, Lz), boundary; min_r = 0.5, temp = 1.0)

coords = [p_info.position for p_info in info.particle_info]
charge = [atoms[p_info.id].charge for p_info in info.particle_info]

γ = (0.5, 0.5)
N_imgs = [1:5:46...]

E_icm = Vector{Vector{Float64}}()
ss = [2.0, 4.0, 6.0]
for s in ss
    t = Vector{Float64}()
    for N_img in N_imgs
        ICMEwald2D_interaction = IcmEwald2DInteraction(n_atoms, s, 2.0, γ, L, N_img)
        energy_ewald = ICMEwald2D_energy(ICMEwald2D_interaction, coords, charge)
        @show N_img, s, energy_ewald
        push!(t, energy_ewald)
    end
    push!(E_icm, t)
end

begin
    E_exact = E_icm[end][end]
    f = Figure(backgroundcolor = RGBf(1.0, 1.0, 1.0), size = (500, 400), fontsize = 20)
    ga = f[1, 1] = GridLayout()
    axl = Axis(ga[1, 1], xlabel = L"N_{img}", ylabel = L"\mathcal{E}_r", yscale = log10)
    for i in 1:length(E_icm)
        lines!(axl, N_imgs[1:end - 1], abs.((E_icm[i][1:end - 1] .- E_exact) ./ E_exact), label = L"s = %$(ss[i])")
        scatter!(axl, N_imgs[1:end - 1], abs.((E_icm[i][1:end - 1] .- E_exact) ./ E_exact))
    end

    axislegend(axl, position = :lb, labelsize = 15)
    save("figs/error_ICM_Ewald2D.png", f)
    save("figs/error_ICM_Ewald2D.pdf", f)
    f
end