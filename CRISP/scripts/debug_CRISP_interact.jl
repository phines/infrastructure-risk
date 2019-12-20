using Glob
using CSV
using DataFrames
using Random
rng = MersenneTwister(100+872)
fold = "casc2/"
l = length(fold);
Files = glob("VACC/results/experiments/mh/"*fold*"res_*n-1.csv")
include("../src/CRISP_Rdist.jl")
events = "data/outage_data/"*fold*"out_case73_noPWS_lx2_n-1"
for path in Files
	## folder of case data
	case = "data/saved_ps/"*path[(33+l):end-4]
	#time steps
	dt = 60 #minutes
	for m in 872
		#save restoration data to folder within results folder:
		filename = "res_out_"*path[(34+l):end-4];
		# run to save csv of resilience cost distribution to the specified out_folder
		#Base line
		fold2 = "casc2_base/"
		out = "/experiments/"*fold2*path[(33+l):end-4]
		out_folder = out*"/$filename-$m.csv"
		if isdir("results"*out)
		else
		    mkdir("results"*out)
		end
		rng = MersenneTwister(100+m);
		comm = false
		nucp = false
		ngi = false
		crt = false
		res = Rdist_interact(m,case,out_folder,events,dt,comm,nucp,ngi,crt)
		#Communications
		fold2 = "casc2+comms/"
		out = "/experiments/"*fold2*path[(33+l):end-4]
		out_folder = out*"/$filename-$m.csv"
		if isdir("results"*out)
		else
		    mkdir("results"*out)
		end
		rng = MersenneTwister(100+m);
		comm = true
		nucp = false
		ngi = false
		crt = false
		res = Rdist_interact(m,case,out_folder,events,dt,comm,nucp,ngi,crt)
		#Nuclear Poissoning
		fold2 = "casc2+nucp/"
		out = "/experiments/"*fold2*path[(33+l):end-4]
		out_folder = out*"/$filename-$m.csv"
		if isdir("results"*out)
		else
		    mkdir("results"*out)
		end
		rng = MersenneTwister(100+m);
		comm = false
		nucp = true
		ngi = false
		crt = false
		res = Rdist_interact(m,case,out_folder,events,dt,comm,nucp,ngi,crt)
		#Natural gas
		fold2 = "casc2+ngi/"
		out = "/experiments/"*fold2*path[(33+l):end-4]
		out_folder = out*"/$filename-$m.csv"
		if isdir("results"*out)
		else
		    mkdir("results"*out)
		end
		rng = MersenneTwister(100+m);
		comm = false
		nucp = false
		ngi = true
		crt = false
		res = Rdist_interact(m,case,out_folder,events,dt,comm,nucp,ngi,crt)
		# Compounding risk over time
		fold2 = "casc2+compi/"
		out = "/experiments/"*fold2*path[(33+l):end-4]
		out_folder = out*"/$filename-$m.csv"
		if isdir("results"*out)
		else
		    mkdir("results"*out)
		end
		rng = MersenneTwister(100+m);
		comm = false
		nucp = false
		ngi = false
		crt = true
		res = Rdist_interact(m,case,out_folder,events,dt,comm,nucp,ngi,crt)
#=		N=m;ps_folder=case;param_file = "";
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
        subgraph = find_subgraphs(ps);
        M = Int64(findmax(subgraph)[1]);
        ps_islands = build_islands(subgraph,ps)
		for i in 1:M
            psi = ps_subset(ps,ps_islands[i]);
            crisp_dcpf_g1_s!(psi);
            # run lsopf
            dt = 10; # first minute, affects ramp rate limits - rateB
            #crisp_lsopf_g1_s!(psi,dt);
            ps.gen.Pg[ps_islands[i].gen] = psi.gen.Pg
            #ps.storage.Ps[ps_islands[i].storage] = psi.storage.Ps
            #ps.storage.E[ps_islands[i].storage] = psi.storage.E
	    	for sh in 1:length(psi.shunt.status)
				if psi.shunt.status[sh] > 1
					psi.shunt.status[sh] = 1.0
				end
	    	end
		end
		#println(case)
		#println(sum(ps.gen.status.==0))
		#println(sum(ps.gen.Pg.==0))
		#println("number of islands")
		#println(M)
		dt = 60
        ti = 60*48
        t0 = 10

		Restore = crisp_Restoration_inter(ps,l_recovery_times,g_recovery_times,dt,
	              ti,t0,gen_on,comm,nucp,ngi,crt)
		outnow = (out_folder[1:end-4])
		CSV.write("results"*outnow*"_restore_nolsopf.csv", Restore)
=#
	end
end
