# Debugging results/case39/*case39_9.*
using Random
rng = MersenneTwister(100);
using CSV
#=
include("src\\CRISP_initiate.jl")
include("src\\CRISP_LSOPF.jl")
include("src\\CRISP_RLSOPF.jl")
include("src\\CRISP_RT.jl")
include("src\\CRISP_network.jl")
=#
include("..\\src\\CRISP_initiate.jl")
include("..\\src\\CRISP_LSOPF.jl")
include("..\\src\\CRISP_RLSOPF.jl")
#include("..\\src\\CRISP_RT.jl")
include("..\\src\\CRISP_network.jl")
## load the case data
for ex in 1:2
ps = import_ps("data/saved_ps/case2736sp_relaxedQ_ps")
ps.branch.status = ones(length(ps.branch.status))
ps.gen.status = ones(length(ps.gen.status))
ps1 = import_ps("data/saved_ps/2736sp_relaxedQ_ps_ex$ex")
ps.shunt = ps.shunt[ps.shunt.P .!=0.0,:]
for s in 1:length(ps.shunt.bus)
    if (sum(ps1.shunt.bus .== ps.shunt.bus[s]) == 0)
        ps.shunt.status[s] = 0
    else
        ps.shunt.status[s] = (ps1.shunt.P[ps1.shunt.bus .== ps.shunt.bus[s]][1]) / ps.shunt.P[s]
    end
end
#crisp_dcpf_g1_s!(ps)
total = sum(ps.shunt.P);
Pd_max = deepcopy(ps.shunt.P);
gen_on = ps.gen.Pg .!= 0;
ps.gen.Pg = ps1.gen.Pg
# parameters
param_file = ""
if isempty(param_file)
    # parameters of distributions for recovery times
    mu_line = 3.66;#lines_dist[2];
    sigma_line = 2.43;#lines_dist[3];
end

#save restoration picture to:
filename = "Polish_Case";
#pull failures from csv file
Lines_Init_State = CSV.File("data/example_outages/lines_tripped_ex$ex.csv") |> DataFrame
lines_state = Lines_Init_State.lines_tripped
l = lines_state[1:end-1]
k = lines_state[2:end]
h = l .== k
h = [false;h]
lines_state = lines_state[.!h]
nl = size(ps.branch,1)
NLines = length(lines_state);
ps.branch.status[Int64.(lines_state)] .= 0;
RecovTimeL = RecoveryTimes(mu_line, sigma_line, NLines);
l_recovery_times = RecTime(RecovTimeL, ps.branch.status).recovery_time;
Gens_Init_State = CSV.File("data/example_outages/gens_tripped_ex$ex.csv") |> DataFrame
ng = size(ps.gen,1);
NGens = length(Gens_Init_State.gens_tripped);
ps.gen.status[Int64.(Gens_Init_State.gens_tripped)] = zeros(NGens);
#RecovTimeG = RecoveryTimes(mu_line, sigma_line, NGens);
g_recovery_times = RecTime(RecovTimeG, ps.gen.status).recovery_time;
g_recovery_times[g_recovery_times .!= 0] = ones(sum(g_recovery_times .!= 0))
## run step 3
dt = 60
ti = 60*48;
t0 = 0
Restore = crisp_RLOPF_v1(ps,l_recovery_times,g_recovery_times,dt,ti,t0,gen_on)
CSV.write("results\\Restoration$(filename)_ex_$ex.csv",Restore)
## run step 4
#ResilienceTri = crisp_res(Restore);
end
