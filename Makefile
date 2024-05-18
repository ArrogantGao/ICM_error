default: init update

init:
	julia --project=. -e 'using Pkg; Pkg.add(url = "https://github.com/HPMolSim/EwaldSummations.jl", branch = "dev");  Pkg.instantiate(); Pkg.precompile();'

update:
	julia --project=. -e 'using Pkg; Pkg.update(); Pkg.precompile();'
