using Random
rng = MersenneTwister(1);
include("../src/CRISP_network.jl");
ps = import_ps("../../../NE_NY_data")
#case = "data/saved_ps/case39_n-1_gen";
folder = "data/saved_ps/black_start";
nodes = ps.bus.id
deg = find_node_degree(nodes,ps.branch.f,ps.branch.t)
reduce_network!(ps,deg)
export_ps(ps,out)
