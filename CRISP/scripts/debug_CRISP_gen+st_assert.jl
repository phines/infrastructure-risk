# test to see why the reference governor code isn't working properly.
#=include("src\\CRISP_initiate.jl")
include("src\\CRISP_LSOPF_gen.jl")
include("src\\CRISP_RLSOPF_gen.jl")
include("src\\CRISP_RT.jl")
include("src\\CRISP_network_gen.jl")=#

include("..\\src\\CRISP_initiate.jl")
include("..\\src\\CRISP_LSOPF_gen.jl")
include("..\\src\\CRISP_RLSOPF_gen.jl")
include("..\\src\\CRISP_RT.jl")
include("..\\src\\CRISP_network_gen.jl")
case1 = "data\\saved_ps\\case73_noPWS+S5\\"
ps = import_ps(case1)
crisp_dcpf_g!(ps)
total = sum(ps.shunt[:P]);
Pd_max = deepcopy(ps.shunt[:P]);
ps0 = deepcopy(ps);

Lines_Init_State = CSV.read("results\\experiments_gen\\7\\res_out_se73_noPWS+S IC48 lines.csv", allowmissing=:none)
Gens_Init_State = CSV.read("results\\experiments_gen\\7\\res_out_se73_noPWS+S IC48 gens.csv", allowmissing=:none)

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
@assert 10^(-6)>=abs(sum(ps.shunt.P)-sum(ps.storage.Ps)-sum(ps.gen.Pg))
for i in 1:M
    psi = ps_subset(ps,ps_islands[i]);
    # run the dcpf
    crisp_dcpf_g_s!(psi);
    # run lsopf
    dt = 10;
    crisp_rlopf_g_s!(psi,Pd_max[ps_islands[i].shunt],dt);
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
@assert 10^(-6)>=abs(sum(ps.shunt.P)-sum(ps.storage.Ps)-sum(ps.gen.Pg))
#Restore = RLSOPF_g_s!(ps,dt,l_failures,g_failures,l_recovery_times,g_recovery_times,Pd_max)
t0 = 10; load_cost=0;
if sum(names(ps.gen).==:time_on) == 0
    ps.gen.time_on = zeros(length(ps.gen.Pg));
end
if sum(names(ps.gen).==:time_off) == 0
    ps.gen.time_off = zeros(length(ps.gen.Pg));
end
if sum(names(ps.gen).==:off) == 0
    ps.gen.off = zeros(length(ps.gen.Pg));
    ps.gen.off[ps.gen.Pg .== 0] .= 1;
end
ps.gen.time_on[ps.gen.off .== 0] .= 5000;
tolerance = 1e-6
if sum(load_cost)==0
    load_cost = ones(length(ps.shunt.P));
end
l_rec_t = l_recovery_times[l_recovery_times.!=0];
g_rec_t = g_recovery_times[g_failures.==0];
g_rec = zeros(length(g_failures));
g_rec[g_failures.==0] = g_rec_t;
Time = sort([l_rec_t; g_rec_t])
dt = 10;
T = 0:dt:maximum(Time)+1000;
total_i = length(T);
load_shed = zeros(total_i+2);
load_shed[1] = 0;
lines_out = zeros(total_i+2);#zeros(length(l_times)+2);
lines_out[2] = length(l_failures) - sum(l_failures);
gens_out = zeros(total_i+2);#zeros(length(l_times)+2);
gens_out[2] = length(g_failures) - sum(g_failures);
# set load shed for the step just before restoration process
load_shed[2] = sum(load_cost.*(Pd_max - ps.shunt[:P]));
time_gen_off = ones(length(ps.gen.Pg));
time_gen_off[ps.gen.Pg .!=0] .=0;
i=1#for i = 1:length(T)
    time = T[i];
    # set failed branches to status 0
    l_failures[time.>=l_recovery_times] .= 1;
    lines_out[i+2] = length(l_failures) - sum(l_failures);
    # apply to network
    ps.branch[:,:status] = l_failures;
    # set newly available generators to status 1;
    current_gen_out = ps.gen.status;
    # update generators operational
    g_failures[time.>=g_rec] .= 1;
    ps.gen.status = g_failures;
    gens_out[i+2] = length(g_failures) - sum(g_failures);
    # update time on and time off
    off = ps.gen.off .!= 0;
    on = ps.gen.off .== 0;
    ps.gen.time_off[off] .+= dt;
    ps.gen.time_on[on] .+= dt;
    #check for islands
    subgraph = find_subgraphs(ps);# add Int64 here hide info here
    M = Int64(findmax(subgraph)[1]);
    if M==1
        crisp_dcpf_g!(ps);
        # run lsopf
        crisp_rlopf_g!(ps,Pd_max);
        crisp_dcpf_g!(ps);
    else
        ps_islands = build_islands(subgraph,ps);# at some point check for changes in islands and don't run power flows if no change
        ## for every island that changed (eventually)

        #for j in 1:M
        j=1;
            psi = ps_subset(ps,ps_islands[j]);
            # run the dcpf
            crisp_dcpf_g!(psi);
            ps.gen[ps_islands[j].gen,:Pg] = psi.gen.Pg
            ps.shunt[ps_islands[j].shunt,:P] = psi.shunt.P
            ps.storage[ps_islands[j].storage,:E] = psi.storage.E
            ps.storage[ps_islands[j].storage,:Ps] = psi.storage.Ps
            ps.branch[ps_islands[j].branch,:Pf] = psi.branch.Pf
            ps.branch[ps_islands[j].branch,:Pt] = psi.branch.Pt
            ps.branch[ps_islands[j].branch,:Qf] = psi.branch.Qf
            ps.branch[ps_islands[j].branch,:Qt] = psi.branch.Qt
            # run lsopf
            crisp_rlopf_g!(psi,Pd_max[ps_islands[j].shunt]);
            # apply the results
            crisp_dcpf_g!(psi);
            ps.gen[ps_islands[j].gen,:Pg] = psi.gen.Pg
            ps.shunt[ps_islands[j].shunt,:P] = psi.shunt.P
            ps.storage[ps_islands[j].storage,:E] = psi.storage.E
            ps.storage[ps_islands[j].storage,:Ps] = psi.storage.Ps
            ps.branch[ps_islands[j].branch,:Pf] = psi.branch.Pf
            ps.branch[ps_islands[j].branch,:Pt] = psi.branch.Pt
            ps.branch[ps_islands[j].branch,:Qf] = psi.branch.Qf
            ps.branch[ps_islands[j].branch,:Qt] = psi.branch.Qt
        #end
    end
    @assert abs(sum(ps.gen.Pg)+sum(ps.storage.Ps)-sum(ps.shunt.P))<=tolerance
    # set load shed for this time step
    load_shed[i+2] = sum(load_cost.*(Pd_max - ps.shunt[:P]));
#end
times = [0.0;t0*1.0;T.+t0*1.0];
perc_load_served = (sum(load_cost.*Pd_max) .- load_shed)./sum(load_cost.*Pd_max);
Restore = DataFrame(time = times, load_shed = load_shed, perc_load_served = perc_load_served, num_lines_out = lines_out, num_gens_out = gens_out);
