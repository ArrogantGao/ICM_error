default: init update

init:
	julia --project=. -e 'using Pkg; Pkg.add(url = "https://github.com/HPMolSim/EwaldSummations.jl");  Pkg.instantiate(); Pkg.precompile();'

update:
	julia --project=. -e 'using Pkg; Pkg.update(); Pkg.precompile();'

run_ICM_Ewald3D:
	touch data/ICM_Ewald3D_2-1.csv
	rm data/ICM_Ewald3D_2-1.csv
	touch data/ICM_Ewald3D_2-1.csv
	echo "N_img,N_slab,γ1,γ2,E_exact,E_ewald,error_r" >> data/ICM_Ewald3D_2-1.csv
	julia --project=. accuracy/ICM_Ewald3D.jl
	julia --project=. drawer/fig_ICM_Ewald3D_2-1.jl