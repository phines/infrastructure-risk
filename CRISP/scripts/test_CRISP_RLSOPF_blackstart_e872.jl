using Glob
using CSV
using DataFrames
using Random
rng = MersenneTwister(1)
fold = "casc2/"
out_fold = "/bs/"
if isdir("results"*out_fold) else mkdir("results"*out_fold) end
fold2 = "casc/"
l = length(fold);
Files = glob("VACC/results/experiments/mh/"*fold*"res_*")
#include("src/CRISP_Rdist_PSCC_comms.jl")
#include("src/CRISP_Rdist_PSCC_casc.jl")
include("../src/CRISP_Rdist.jl")
include("../src/CRISP_RLSOPF.jl")
include("../src/CRISP_network.jl")
#events = "data/outage_data/out_case73_noPWS_lx2_n-1"
#events = "data/outage_data/communication_factor/out_case73_noPWS_lx2_n-1"

events = "data/outage_data/"*fold*"out_case73_noPWS_lx2_n-1"
comm = false;
nucp = false;
ngi = false;
crt = false;
dt = 60*10;
N = 872;
Lines_Init_State = CSV.File(events*"_lines$N.csv") |> DataFrame
Gens_Init_State = CSV.File(events*"_gens$N.csv") |> DataFrame
ps_folder = "data/saved_ps/bs/case73_noPWS_lx2_n-1"
out_folder = out_fold*"res_case73_noPWS_lx2_n-1.csv"
ps = import_ps("$ps_folder")
ng = size(ps.gen,1);
ps.gen[!,:state] = Vector{Enum}(undef,ng);
ps.gen[!,:time_in_state] = zeros(length(ps.gen.bus));
set_gen_states!(ps);
crisp_dcpf_g1_s!(ps)
total = sum(ps.shunt.P);
Pd_max = deepcopy(ps.shunt.P);
gen_on = ps.gen.Pg .!= 0;
l_failures = Lines_Init_State.state;
ps.branch.status[l_failures .== 0] .= 0;
l_recovery_times = deepcopy(Lines_Init_State.recovery_time);
# generator states and recovery times
g_failures = ones(size(ps.gen,1));
g_failures[1:length(Gens_Init_State.state)] = Gens_Init_State.state;
ps.gen.status[g_failures .== 0] .= 0;
ps.gen.Pg[g_failures .== 0] .= 0;
ps.gen.state[g_failures .== 0] .= Damaged;
ps.gen.time_in_state[g_failures .== 0] .= 0;
gens_recovery_time = zeros(size(ps.gen,1));
gens_recovery_time[1:length(Gens_Init_State.recovery_time)] = Gens_Init_State.recovery_time;
t = 60
subgraph = find_subgraphs(ps);# add Int64 here hide info here
M = Int64(findmax(subgraph)[1]);
ps_islands = build_islands(subgraph,ps)
Time = 0:60:1000
t_window = dt
ul = line_stats(ps,Time,t_window,l_recovery_times)
# varying load over the course of the optimization
Pg_max = vary_gen_cap(ps,Time,t_window)
Pd_max = vary_load(ps,Time,t_window)
load_cost = ones(length(ps.shunt.bus))
for j in 1:M
    psi = ps_subset(ps,ps_islands[j])
    i_subset = 1:1+1
    uli = ul[ps_islands[j].branch,i_subset]
    Pd_maxi = Pd_max[ps_islands[j].shunt,i_subset]
    Pg_maxi = Pg_max[ps_islands[j].gen,i_subset]
    crisp_lsopf_bs!(psi,dt,uli,Pd_maxi,Pg_maxi,load_cost[ps_islands[j].shunt])
    ps.gen.Pg[ps_islands[j].gen] = psi.gen.Pg
    ps.storage.Ps[ps_islands[j].storage] = psi.storage.Ps
    ps.storage.E[ps_islands[j].storage] = psi.storage.E
    ps.shunt.status[ps_islands[j].shunt] = psi.shunt.status
end
println(ps.shunt.status)
println(ps.gen.state)
println(gens_recovery_time)
pgmax,ps,gens_recovery_time = black_start_gen_cap!(ps,t,dt,gens_recovery_time)
println(ps.gen.state)
println(gens_recovery_time)
#=
for path in Files
	## folder of case data
	case = "data/saved_ps/"*path[(33+l):end-4]
	out = "/experiments/"*fold*path[(33+l):end-4]
	#time steps
	dt = 60 #minutes
	for m in 872
		#save restoration data to folder within results folder:
		filename = "res_out_"*path[(34+l):end-4];
		out_folder = out*"/$filename-$m.csv"
		if isdir("results"*out)
		else
		    mkdir("results"*out)
		end
		# run to save csv of resilience cost distribution to the specified out_folder
		#res = Resilience(m,case,out_folder,events,dt)
		N=m;ps_folder=case;param_file = "";

		#constants
	    debug=1;
	    tolerance1 = 10^(-4);
	    Num = 1; #number of events to run
	    ## N = which failure scenario
	    # initialize vector of costs from events
	    NumLinesOut = Array{Float64}(undef,Num,1);
	    LoadShed0 =  Array{Float64}(undef,Num,1);
	    MaxRestorationTime = Array{Float64}(undef,Num,1);
	    LoadServedTime = Array{Float64}(undef,Num,1);
	    ResilienceTri = Array{Float64}(undef,Num,1);
	    ## load the case data
	    ps = import_ps("$ps_folder")
	    ps.shunt = ps.shunt[ps.shunt.P .!=0.0,:]
	    crisp_dcpf_g1_s!(ps)
	    total = sum(ps.shunt.P);
	    Pd_max = deepcopy(ps.shunt.P);
	    gen_on = ps.gen.Pg .!= 0;
	    ps0 = deepcopy(ps)
		iterat = 1
		# step 1
        Lines_Init_State = CSV.File(events*"_lines$N.csv") |> DataFrame
        Gens_Init_State = CSV.File(events*"_gens$N.csv") |> DataFrame
        l_failures = Lines_Init_State.state;
        ps.branch.status[l_failures .== 0] .= 0;
        NumLinesOut[iterat] = length(l_failures) - sum(l_failures)
        l_recovery_times = Lines_Init_State.recovery_time;
        # generator states and recovery times
        g_failures = ones(size(ps.gen,1));
        g_failures[1:length(Gens_Init_State.state)] = Gens_Init_State.state;
        ps.gen.status[g_failures .== 0] .= 0;
        ps.gen.Pg[g_failures .== 0] .= 0;
        g_recovery_times = zeros(size(ps.gen,1));
        g_recovery_times[1:length(Gens_Init_State.recovery_time)] = Gens_Init_State.recovery_time;
		#check for islands
        #subgraph = find_subgraphs(ps);
        #M = Int64(findmax(subgraph)[1]);
        #ps_islands = build_islands(subgraph,ps)
		#for i in 1:M
            #psi = ps_subset(ps,ps_islands[i]);
            #crisp_dcpf_g1_s!(psi);
            # run lsopf
            #dt = 10; # first minute, affects ramp rate limits - rateB
            #crisp_lsopf_g1_s!(psi,dt);
            #ps.gen.Pg[ps_islands[i].gen] = psi.gen.Pg
            #ps.storage.Ps[ps_islands[i].storage] = psi.storage.Ps
            #ps.storage.E[ps_islands[i].storage] = psi.storage.E
	    	#for sh in 1:length(psi.shunt.status)
				#if psi.shunt.status[sh] > 1
					#psi.shunt.status[sh] = 1.0
				#end
	    	#end
		#end
		#println(case)
		#println(sum(ps.gen.status.==0))
		#println(sum(ps.gen.Pg.==0))
		#println("number of islands")
		#println(M)
		dt = 60
        ti = 60*48;
        t0 = 10
		Restore = crisp_Restoration_var(ps,l_recovery_times,g_recovery_times,dt,ti,t0,gen_on)
		outnow = (out_folder[1:end-4])
		CSV.write("results"*outnow*"_restore_newlsopf_obj_Pg.csv", Restore)
	end
end
=#
