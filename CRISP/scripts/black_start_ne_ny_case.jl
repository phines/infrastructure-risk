using Random
rng = MersenneTwister(1);
include("../src/CRISP_network.jl");
case = "../../../NE_NY_data/ps_original"
#case = "data/saved_ps/case39_n-1_gen";
folder = "../../../NE_NY_data/ps_reduced_bs"
if isdir(folder) else mkdir(folder) end
out = folder
ps0 = import_ps(case);
ps = deepcopy(ps0);
ps.shunt = ps.shunt[ps.shunt.P .!=0.0,:];
ng = size(ps.gen,1);
ps.gen[!,:service_load] = 0.05.*ps.gen.Pmax;
ps.gen[!,:sl_status] = ones(ng);
choose_gens_black_start!(ps,0.5, 20);
ps.gen.Pmax[.!ps.gen.black_start] .+=0.05.*ps.gen.Pmax[.!ps.gen.black_start];
export_ps(ps,out)
