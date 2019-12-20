# Debugging results/case39/*case39_9.*
rng = MersenneTwister(100);
using CSV
include("..\\src\\CRISP_initiate.jl")
include("..\\src\\CRISP_LSOPF.jl")
include("..\\src\\CRISP_RLSOPF.jl")
include("..\\src\\CRISP_RT.jl")
include("..\\src\\CRISP_network.jl")
## load the case data
ps = import_ps("data/saved_ps/case2736sp_relaxedQ_ps")
crisp_dcpf!(ps)
total = sum(ps.shunt[:P]);
Pd_max = deepcopy(ps.shunt[:P]);

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
Lines_Init_State = CSV.File("data/example_outages/lines_tripped_ex1.csv") |> DataFrame
nl = size(ps.branch,1)
NLines = length(Lines_Init_State.lines_tripped);
ps.branch.status[Int64.(Lines_Init_State.lines_tripped)] = zeros(NLines);
RecovTimeL = RecoveryTimes(mu_line, sigma_line, NLines);
l_recovery_times = RecTime(RecovTimeL, ps.branch.status);
Gens_Init_State = CSV.File("data/example_outages/gens_tripped_ex1.csv") |> DataFrame
ng = size(ps.gen,1);
NGens = length(Gens_Init_State.gens_tripped);
ps.gen.status[Int64.(Gens_Init_State.gens_tripped)] = zeros(NGens);
RecovTimeG = RecoveryTimes(mu_line, sigma_line, NGens);
g_recovery_times = RecTime(RecovTimeG, ps.gen.status);
## run step 3
dt = 60
ti = 60*48;
t0 = 10
Restore = crisp_RLOPF_v1(ps,l_recovery_times,g_recovery_times,dt,t_window,t0,gen_on)
CSV.write("results\\Restoration$filename.csv",Restore)
## run step 4
ResilienceTri = crisp_res(Restore);
