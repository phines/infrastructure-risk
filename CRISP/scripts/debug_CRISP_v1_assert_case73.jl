# test to see why the code isn't working properly.
#=include("src\\CRISP_initiate.jl")
include("src\\CRISP_LSOPF1.jl")
include("src\\CRISP_RLSOPF.jl")
include("src\\CRISP_RT.jl")
include("src\\CRISP_network.jl")=#

include("..\\src\\CRISP_initiate.jl")
include("..\\src\\CRISP_LSOPF1.jl")
include("..\\src\\CRISP_RLSOPF.jl")
include("..\\src\\CRISP_RT.jl")
include("..\\src\\CRISP_network.jl")

ps = import_ps("data\\saved_ps\\case73_noPWS_n-1\\")
crisp_dcpf!(ps)
total = sum(ps.shunt[:P]);
Pd_max = deepcopy(ps.shunt[:P]);
ps0 = deepcopy(ps);

Lines_Init_State = CSV.read("results\\100\\case73\\res_out_case39_n-1_p1 IC544.csv", allowmissing=:none)

state = Lines_Init_State[:,1];
ps.branch.status[state .==0] .= 0;
NumLinesOut = length(state) - sum(state)
recovery_times = Lines_Init_State[:,2];
MaxRestorationTime = maximum(recovery_times);
failures = state;
#check for islands
subgraph = find_subgraphs(ps);
M = Int64(findmax(subgraph)[1]);
ps_islands = build_islands(subgraph,ps)

for i in 1:M
    println(i)
    psi = ps_subset(ps,ps_islands[i]);
    # run the dcpf
    crisp_dcpf!(psi);
    # run lsopf
    crisp_lsopf1!(psi);
    ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
    ps.shunt[ps_islands[i].shunt,:P] = psi.shunt.P
    crisp_dcpf!(psi);
end
@assert 10^(-6)>=abs(sum(ps.shunt.P)-sum(ps.gen.Pg))
Restore = RLSOPF!(total,ps,state,recovery_times,Pd_max)
