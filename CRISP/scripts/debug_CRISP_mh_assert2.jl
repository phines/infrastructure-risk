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
case1 = "data\\saved_ps\\case73_noPWS+S5\\"
ps = import_ps(case1)
ps.shunt = ps.shunt[ps.shunt.P .!=0.0,:]
crisp_dcpf_g_s!(ps)
total = sum(ps.shunt[:P]);
Pd_max = deepcopy(ps.shunt[:P]);
#add columns to keep track of the time each generator is on or off
if sum(names(ps.gen).==:time_on) == 0
    ps.gen.time_on = zeros(length(ps.gen.Pg));
end
if sum(names(ps.gen).==:time_off) == 0
    ps.gen.time_off = zeros(length(ps.gen.Pg));
    ps.gen.time_off[ps.gen.Pg.==0] .= ps.gen.minDownTimeHr[ps.gen.Pg.==0];
end
ps0 = deepcopy(ps);

Lines_Init_State = CSV.File("results\\experiments_gen_stor\\case73\\mh\\res_out_case73_noPWS+S5 IC1 lines.csv")
Gens_Init_State = CSV.File("results\\experiments_gen_stor\\case73\\mh\\res_out_case73_noPWS+S5 IC1 gens.csv")

l_failures = Lines_Init_State.state;
ps.branch.status[l_failures .==0] .= 0;
NumLinesOut = length(l_failures) - sum(l_failures)
l_recovery_times = Lines_Init_State.recovery_time;
# generator states and recovery times
g_failures = Gens_Init_State.state;
ps.gen.status[g_failures.==0] .= 0;
ps.gen.Pg[g_failures.==0] .= 0;
g_recovery_times = Gens_Init_State.recovery_time;
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
Restore = crisp_Restore_mh(ps,l_recovery_times,g_recovery_times,dt,t_window,t0)
