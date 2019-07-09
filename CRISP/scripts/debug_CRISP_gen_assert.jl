# test to see why the reference governor code isn't working properly.
include("..\\src\\CRISP_initiate.jl")
include("..\\src\\CRISP_LSOPF_gen.jl")
include("..\\src\\CRISP_RLSOPF_gen.jl")
include("..\\src\\CRISP_RT.jl")
include("..\\src\\CRISP_network_gen.jl")

ps = import_ps("data\\case73")
crisp_dcpf_g!(ps)
total = sum(ps.shunt[:P]);
Pd_max = deepcopy(ps.shunt[:P]);
ps0 = deepcopy(ps);

Lines_Init_State = CSV.read("results\\experiments_gen\\3\\res_out_case73_gendata IC21 lines.csv", allowmissing=:none)
Gens_Init_State = CSV.read("results\\experiments_gen\\3\\res_out_case73_gendata IC21 gens.csv", allowmissing=:none)

state = Lines_Init_State[:,1];
ps.branch.status[state .==0] .= 0;
NumLinesOut = length(state) - sum(state)
recovery_times = Lines_Init_State[:,2];
MaxRestorationTime = maximum(recovery_times);
failures = state;
# generator states and recovery times
gens_state = Gens_Init_State[:,1];
ps.gen.status[gens_state.==0] .= 0;
ps.gen.Pg[gens_state.==0] .= 0;
gens_recovery_time = Gens_Init_State[:,2];
#check for islands
subgraph = find_subgraphs(ps);
M = Int64(findmax(subgraph)[1]);
ps_islands = build_islands(subgraph,ps) #=;

for i in 1:M
    psi = ps_subset(ps,ps_islands[i]);
    # run the dcpf
    crisp_dcpf_g!(psi);
    # run lsopf
    crisp_lsopf_g!(psi);
    ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
    ps.shunt[ps_islands[i].shunt,:P] = psi.shunt.P
    crisp_dcpf_g!(psi);
end
@assert 10^(-6)>=abs(sum(ps.shunt.P)-sum(ps.gen.Pg))
Restore = RLSOPF_g!(ps,state,gens_state,recovery_times,gens_recovery_time,Pd_max)
=#
