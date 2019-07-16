# test to see why the reference governor code isn't working properly.
#=include("src\\CRISP_initiate.jl")
include("src\\CRISP_LSOPF_gen1.jl")
include("src\\CRISP_RLOPF_movh.jl")
include("src\\CRISP_RT.jl")
include("src\\CRISP_network_gen.jl")=#

include("..\\src\\CRISP_initiate.jl")
include("..\\src\\CRISP_LSOPF_gen1.jl")
include("..\\src\\CRISP_RLOPF_movh.jl")
include("..\\src\\CRISP_RT.jl")
include("..\\src\\CRISP_network_gen.jl")
case1 = "data\\saved_ps\\case39_n-1_gen\\"
ps = import_ps(case1)
ps.shunt = ps.shunt[ps.shunt.P .!=0.0,:]
crisp_dcpf_g_s!(ps)
total = sum(ps.shunt[:P]);
Pd_max = deepcopy(ps.shunt[:P]);
ps0 = deepcopy(ps);

Lines_Init_State = CSV.read("results\\experiments_gen_mvh\\res_out_case39_n-1_gen IC3 lines.csv", allowmissing=:none)
Gens_Init_State = CSV.read("results\\experiments_gen_mvh\\res_out_case39_n-1_gen IC3 gens.csv", allowmissing=:none)

l_failures = Lines_Init_State[:,1];
ps.branch.status[l_failures .==0] .= 0;
NumLinesOut = length(l_failures) - sum(l_failures)
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
dt = 1;
for i in 1:M
    psi = ps_subset(ps,ps_islands[i]);
    crisp_lsopf_g_s!(psi,dt);
    ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
    ps.storage[ps_islands[i].storage,:Ps] = psi.storage.Ps
    ps.storage[ps_islands[i].storage,:E] = psi.storage.E
    ps.shunt[ps_islands[i].shunt,:status] = psi.shunt.status
    println("Pg = ")
    println(sum(psi.gen.Pg))
    println("P = ")
    println(sum(psi.shunt.P .* psi.shunt.status))
    println("Ps = ")
    println(sum(psi.storage.Ps))
    @assert 10^(-4)>=abs(sum(psi.shunt.P .* psi.shunt.status)-sum(psi.storage.Ps)-sum(psi.gen.Pg))
    #=
    # run the dcpf
    crisp_dcpf_g_s!(psi);
    # run lsopf
    dt = 10;
    crisp_lsopf_g_s1!(psi,dt);
    ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
    ps.storage[ps_islands[i].storage,:Ps] = psi.storage.Ps
    ps.storage[ps_islands[i].storage,:E] = psi.storage.E
    ps.shunt[ps_islands[i].shunt,:P] = psi.shunt.P
    crisp_dcpf_g_s!(psi);
    ps.branch[ps_islands[i].branch,:Pf] = psi.branch.Pf
    ps.branch[ps_islands[i].branch,:Pt] = psi.branch.Pt
    ps.branch[ps_islands[i].branch,:Qf] = psi.branch.Qf
    ps.branch[ps_islands[i].branch,:Qt] = psi.branch.Qt
    =#
end
println("Pg = ")
println(sum(ps.gen.Pg))
println("P = ")
println(sum(ps.shunt.P .* ps.shunt.status))
println("Ps = ")
println(sum(ps.storage.Ps))
@assert 10^(-4)>=abs(sum(ps.shunt.P .* ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg))

#LoadShed0[iterat] = total-sum(ps.shunt.P);
## run step 3
dt = 15;
t_window = dt;#10
t0 = 10
#crisp_mh_rlopf!(ps,dt,time)
Restore = crisp_Restore(ps,l_recovery_times,g_recovery_times,dt,t_window,t0)
