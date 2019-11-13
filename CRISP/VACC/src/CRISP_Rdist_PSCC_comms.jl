using CSV
include("CRISP_LSOPF_gen1.jl")
include("CRISP_RLOPF_mh_varLoad.jl")
include("CRISP_RT.jl")
include("CRISP_network_gen.jl")

function Resilience(N,ps_folder,out_folder,events,dt;param_file = "")
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
    crisp_dcpf_g_s!(ps)
    total = sum(ps.shunt.P);
    Pd_max = deepcopy(ps.shunt.P);
    gen_on = ps.gen.Pg .!= 0;
    ps0 = deepcopy(ps);
    for iterat in 1:Num
        ps = deepcopy(ps0); # reset network to original state
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
            crisp_dcpf_g_s!(psi);
            # run lsopf
            dt = 10; # first minute, affects ramp rate limits - rateB
            crisp_lsopf_g_s!(psi,dt);
            ps.gen.Pg[ps_islands[i].gen] = psi.gen.Pg
            ps.storage.Ps[ps_islands[i].storage] = psi.storage.Ps
            #ps.storage.E[ps_islands[i].storage] = psi.storage.E
	    for sh in 1:length(psi.shunt.status)
		if psi.shunt.status[sh] > 1
			psi.shunt.status[sh] = 1.0
		end
	    end
            ps.shunt.status[ps_islands[i].shunt] = psi.shunt.status
        end
        println(iterat)
        @assert total>=sum(ps.shunt.P .* ps.shunt.status)
        @assert 2M*tolerance1 >= abs(sum(ps.shunt.P .*ps.shunt.status)-sum(ps.gen.Pg)-sum(ps.storage.Ps))
        tolerance = 10^(-10);
        LoadShed0[iterat] = total-sum(ps.shunt.P .* ps.shunt.status);
        ## run step 3
        dt = 60
        ti = 60*48;
        t0 = 10
        Restore = crisp_Restoration(ps,l_recovery_times,g_recovery_times,dt,ti,t0,gen_on)
        if debug==1
            outnow = (out_folder[1:end-4]);
            CSV.write("results"*outnow*"_restore.csv", Restore)
        end
        ## find the time to restore the grid to 99.9% load served
        K = abs.(Restore.perc_load_served .- 1) .<= 0.001;
        K[1] = false;
        if isempty( Restore.time[K]) && (size(Restore)[1] > 1)
            error("Never solved to full restoration.")
        elseif isempty( Restore.time[K]) && (size(Restore)[1] == 1)
            LoadServedTime = t0;
        else
            R = Restore.time[K];
            LoadServedTime[iterat] = R[1] - Restore.time[2];
            ## run step 4
            ResilienceTri[iterat] = crisp_res(Restore);
            println(ResilienceTri)
        end
    end
    ## make dataframes
    case_res = DataFrame(resilience = ResilienceTri[:,1]);
    #case_lines = DataFrame(lines_out = NumLinesOut[:,1], no_load_shed_time = LoadServedTime[:,1], initial_load_shed = LoadShed0[:,1]);#, full_rest_time = MaxRestorationTime[:,1]);

    ## save data
    CSV.write("results/$out_folder", case_res);
    #CSV.write("results/$(out_folder[1:end-4])_lines.csv", case_lines);
end

#=
using CSV
include("VACC\\src\\CRISP_initiate.jl")
include("VACC\\src\\CRISP_LSOPF_gen1.jl")
include("VACC\\src\\CRISP_RLOPF_mh_varLoad.jl")
include("VACC\\src\\CRISP_RT.jl")
include("VACC\\src\\CRISP_network_gen.jl")
=#
