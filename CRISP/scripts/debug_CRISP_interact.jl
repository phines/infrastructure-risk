using Glob
using CSV
using DataFrames
using Random
rng = MersenneTwister(100+872+rand(collect(1:100000000),1)[1])
fold = "casc2/"
l = length(fold);
Files = glob("VACC/results/experiments/mh/"*fold*"res_*n-1.csv")
events = "data/outage_data/"*fold*"out_case73_noPWS_lx2_n-1"
#=
include("src/CRISP_Rdist.jl")
include("src/CRISP_RLSOPF.jl")
include("src/CRISP_interact.jl")
=#
include("../src/CRISP_Rdist.jl")
include("../src/CRISP_RLSOPF.jl")
include("../src/CRISP_interact.jl")
path=Files[1]
#for path in Files
	## folder of case data
	case = "data/saved_ps/"*path[(33+l):end-4]
	#time steps
	dt = 60 #minutes
	m = 872;
	#for m in 872
		#save restoration data to folder within results folder:
		filename = "res_out_"*path[(34+l):end-4];
		#=
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
=#
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
		nucp = true
		ngi = false
		crt = false
#		res = Rdist_interact(m,case,out_folder,events,dt,comm,nucp,ngi,crt)
		N=m;ps_folder=case;param_file = "";
		#constants
	    debug=1;
	    tolerance1 = 10^(-4);
	    Num = 1; #number of events to run
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
        NumLinesOut = length(l_failures) - sum(l_failures)
        l_recovery_times = Lines_Init_State.recovery_time;
        # generator states and recovery times
        g_failures = ones(size(ps.gen,1));
        g_failures[1:length(Gens_Init_State.state)] = Gens_Init_State.state;
        ps.gen.status[g_failures .== 0] .= 0;
        ps.gen.Pg[g_failures .== 0] .= 0;
        g_recovery_times = zeros(size(ps.gen,1));
        g_recovery_times[1:length(Gens_Init_State.recovery_time)] = Gens_Init_State.recovery_time;
		dt = 60
        ti = 60*48
        t0 = 10
		Restore = crisp_RLOPF_inter(ps,l_recovery_times,g_recovery_times,dt,
	              ti,t0,gen_on,comm,nucp,ngi,crt)
		#=t_window = ti;
		load_cost=0;com_bl_a=4;com_bl_b = 24;c_factor=1.5;comp_t=8*60;factor=1.5;
		tolerance = 10^(-6);
	    if sum(load_cost)==0
	        load_cost = ones(length(ps.shunt.P));
	    end
	    ti = t0;
	    if comm
	        comm_battery_limits = comm_battery_lim(size(ps.bus,1),com_bl_a,com_bl_b)
	    end
	    recTime = maximum([maximum(l_recovery_times) maximum(g_recovery_times)]);
	    # set time line
	    EndTime = (t0+recTime+((maximum(ps.gen.minDownTimeHr)+maximum(ps.gen.minUpTimeHr))*60));
	    Time = t0:dt:EndTime
	    #save initial values
	    load_shed = sum(load_cost.*(ps.shunt.P - ps.shunt.P.*ps.shunt.status));
	    perc_load_served = (sum(load_cost.*ps.shunt.P) .- load_shed)./sum(load_cost.*ps.shunt.P);
	    lines_out = length(ps.branch.status) - sum(ps.branch.status);
	    gens_out = length(ps.gen.status) - sum(ps.gen.status);
	    Restore = DataFrame(time = ti, load_shed = load_shed, perc_load_served = perc_load_served,
	    lines_out = lines_out, gens_out = gens_out)
	    cv = deepcopy(Restore);
	    # find generator status
	    ug = gen_on_off(ps,Time,t_window,gen_on,g_recovery_times)
	    # find line status
	    ul = line_stats(ps,Time,t_window,l_recovery_times)
	    # varying load over the course of the optimization
	    Pd_max = vary_load(ps,Time,t_window)
	    # varying generation capacity over the optimization
	    Pg_max = vary_gen_cap(ps,Time,t_window)
		#i=1
		for i in 1:length(Time)
	        if nucp
	            Pg_i = deepcopy(ps.gen.Pg);
	        end
	        # update time
	        ti = Time[i]-t0;
	        # remove failures as the recovery time is reached
	        ps.branch.status[ti .>= l_recovery_times] .= 1;
	        #comm_count[ti .>= l_recovery_times] .= 100;
	        ps.gen.status[ti .>= g_recovery_times] .= 1;
	        # find the number of islands in ps
	        subgraph = find_subgraphs(ps);# add Int64 here hide info here
	        M = Int64(findmax(subgraph)[1]);
	        ps_islands = build_islands(subgraph,ps)
	        for j in 1:M
	            psi = ps_subset(ps,ps_islands[j])
	            i_subset = i:i+1
	            ugi = ug[ps_islands[j].gen,i_subset]
	            uli = ul[ps_islands[j].branch,i_subset]
	            Pd_maxi = Pd_max[ps_islands[j].shunt,i_subset]
	            Pg_maxi = Pg_max[ps_islands[j].gen,i_subset]
	            crisp_mh_lsopf_var!(psi,dt,ugi,uli,Pd_maxi,Pg_maxi,load_cost[ps_islands[j].shunt])
	            ps.gen.Pg[ps_islands[j].gen] = psi.gen.Pg
	            ps.storage.Ps[ps_islands[j].storage] = psi.storage.Ps
	            ps.storage.E[ps_islands[j].storage] = psi.storage.E
	            ps.shunt.status[ps_islands[j].shunt] = psi.shunt.status
	        end
	        if nucp
	            g_recovery_times = nuclear_poissoning(ps,Pg_i,g_recovery_times,ti)
	        end
	        if comm
	            if (ti >= com_bl_a*60) .& (ti<= com_bl_b*60) #most communcation towers have batteries which have a capacity to cover from 4 to 24 hour
	                l_recovery_times = communication_interactions!(ps,l_recovery_times,comm_battery_limits,ti,c_factor)
	            end
	        end
	        if crt
	            if (abs(ti./comp_t - round(ti./comp_t)) <= tolerance) & (ti > 0)
	                l_recovery_times = compound_rest_times!(ps,l_recovery_times,factor,ti)
	            end
	        end
	        # save current values
	        cv.time .= ti+t0;
	        cv.load_shed .= sum(load_cost.*(Pd_max[:,i+1] - Pd_max[:,i+1].*ps.shunt.status));
	        cv.perc_load_served .= (sum(load_cost.*Pd_max[:,i+1]) .- cv.load_shed)./sum(load_cost.*Pd_max[:,i+1]);
	        cv.lines_out .= length(ps.branch.status) - sum(ps.branch.status);
	        cv.gens_out .= length(ps.gen.status) - sum(ps.gen.status);
	        append!(Restore,cv)
	        @assert 10^(-4)>=abs(sum(Pd_max[:,i+1] .* ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg))
	    end =#
		outnow = (out_folder[1:end-4])
		CSV.write("results"*outnow*"_restore_nucp.csv", Restore)
	#end
#end
