using Random
rng = MersenneTwister(1);
include("../src/CRISP_network.jl");
case = "data/saved_ps/case73_noPWS_lx2_n-1";
#case = "data/saved_ps/case39_n-1_gen";
folder = "data/saved_ps/black_start";
if isdir(folder) else mkdir(folder) end
out = folder*case[14:end]
ps0 = import_ps(case);
ps = deepcopy(ps0);
ps.shunt = ps.shunt[ps.shunt.P .!=0.0,:];
ng = size(ps.gen,1);
@enum GenState OutOfOpperation Damaged Off ShuttingDown WarmingUp On
ps.gen[!,:state] = Vector{Enum}(undef,ng);
ps.gen[!,:time_in_state] = zeros(length(ps.gen.bus));
ps.gen[!,:service_load] = 0.05.*ps.gen.Pmax;
set_gen_states!(ps);
choose_gens_black_start!(ps,0.5, 20);
out
#export_ps(ps,out)
