using ExTinyMD, EwaldSummations
using CSV, DataFrames
using Random
Random.seed!(1234)

n_atoms = 39
Lxy = 15.0
Lz = 5.0
L = (Lxy, Lxy, Lz)
boundary = ExTinyMD.Q2dBoundary(Lxy, Lxy, Lz)

atoms = Vector{Atom{Float64}}()
for i in 1:26
    push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
end

for i in 1:13
    push!(atoms, Atom(type = 2, mass = 1.0, charge = - 2.0))
end

info = SimulationInfo(n_atoms, atoms, (0.0, Lxy, 0.0, Lxy, 0.5, Lz - 0.5), boundary; min_r = 0.5, temp = 1.0)

coords = [p_info.position for p_info in info.particle_info]
charge = [atoms[p_info.id].charge for p_info in info.particle_info]

Δs = [0.6, 0.95, 1.0]
N_imgs = [2:2:16...]
N_slabs = [3, 5, 7]

data_file = "data/ICM_Ewald3D_2-1.csv"

for Δ in Δs
    for γ in [(Δ, Δ), (-Δ, -Δ), (-Δ, Δ)]
        ICMEwald2D_Interaction = IcmEwald2DInteraction(n_atoms, 6.0, 1.0, γ, L, 20)
        energy_exact = ICMEwald2D_energy(ICMEwald2D_Interaction, coords, charge)
        for n_slab in N_slabs
            for N_img in N_imgs
                N_slab = (3 * n_slab - 1) ÷ 2
                ICMEwald3D_interaction = IcmEwald3DInteraction(n_atoms, 4.0, 1.0, γ, L, N_img, N_slab)
                energy_ewald = ICMEwald3D_energy(ICMEwald3D_interaction, coords, charge)
                er = abs(energy_ewald - energy_exact) / abs(energy_exact)
                @show N_img, N_slab, γ, er
                df = DataFrame(N_img = N_img, N_slab = n_slab, γ1 = γ[1], γ2 = γ[2], E_exact = energy_exact, E_ewald = energy_ewald, error_r = er)
                CSV.write(data_file, df, append = true)
            end
        end
    end
end