using Random
rng = MersenneTwister(1)
include("../src/CRISP_network.jl")
case = "data/saved_ps/case73_noPWS_lx2_n-1"
#case = "data/saved_ps/case39_n-1_gen"
ps0 = import_ps(case)
ps = deepcopy(ps0)
ps.shunt = ps.shunt[ps.shunt.P .!=0.0,:]
ng = size(ps.gen,1)
@enum GenState OutOfOpperation Damaged Off ShuttingDown WarmingUp On
ps.gen[!,:state] = Vector{Enum}(undef,ng)
gen_states!(ps)
choose_gens_black_start!(ps,0.5, 20)
println(minimum(ps.gen.Pmax))

#export_ps(ps)
