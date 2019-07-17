include("..\\src\\CRISP_network_gen.jl")
include("..\\src\\CRISP_initiate.jl")
include("..\\src\\CRISP_LSOPF_gen1.jl")
include("..\\src\\CRISP_RT.jl")
include("..\\src\\CRISP_RLOPF_movh.jl")
case = "data\\saved_ps\\case73_noPWS"
ps = import_ps(case)
ps.shunt = ps.shunt[ps.shunt.P .!=0.0,:]
crisp_dcpf_g_s!(ps)
total = sum(ps.shunt[:P]);
Pd_max = deepcopy(ps.shunt[:P]);
ps0 = deepcopy(ps);


Lines_Init_State = CSV.read("results\\experiments_gen\\7\\res_out_se73_noPWS+S IC48 lines.csv", allowmissing=:none)
Gens_Init_State = CSV.read("results\\experiments_gen\\7\\res_out_se73_noPWS+S IC48 gens.csv", allowmissing=:none)

l_failures = Lines_Init_State[:,1];
ps.branch.status[l_failures .==0] .= 0;
NumLinesOut = length(l_failures) - sum(l_failures);
l_recovery_times = Lines_Init_State[:,2];
# generator states and recovery times
g_failures = Gens_Init_State[:,1];
ps.gen.status[g_failures.==0] .= 0;
ps.gen.Pg[g_failures.==0] .= 0;
g_recovery_times = Gens_Init_State[:,2];
#check for islands
subgraph = find_subgraphs(ps);
M = Int64(findmax(subgraph)[1]);
ps_islands = build_islands(subgraph,ps)
for i in 1:M
    psi = ps_subset(ps,ps_islands[i]);
    # run the dcpf
    crisp_dcpf_g_s!(psi);
    # run lsopf
    dt = 10;
    crisp_lsopf_g_s!(psi,dt);
    ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
    ps.storage[ps_islands[i].storage,:Ps] = psi.storage.Ps
    ps.storage[ps_islands[i].storage,:E] = psi.storage.E
    ps.shunt[ps_islands[i].shunt,:P] = psi.shunt.P
    crisp_dcpf_g_s!(psi);
    ps.branch[ps_islands[i].branch,:Pf] = psi.branch.Pf
    ps.branch[ps_islands[i].branch,:Pt] = psi.branch.Pt
    ps.branch[ps_islands[i].branch,:Qf] = psi.branch.Qf
    ps.branch[ps_islands[i].branch,:Qt] = psi.branch.Qt
end
dt = 1
time = 10
t0 = 10
println(ps.shunt.P)
println(ps.shunt.status)
#crisp_mh_rlopf!(ps,dt,time)
crisp_Restore_mh(ps,l_recovery_times,g_recovery_times,dt,time,t0)
