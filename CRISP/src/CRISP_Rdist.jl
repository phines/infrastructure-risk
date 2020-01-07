using CSV
include("CRISP_initiate.jl")
include("CRISP_LSOPF.jl")
include("CRISP_RLSOPF.jl")
include("CRISP_RT.jl")
include("CRISP_network.jl")

function Rdist_interact(N,ps_folder,out_folder,events,dt,comm,nucp,ngi,crt)
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
    # step 1
    Lines_Init_State = CSV.File(events*"_lines$N.csv") |> DataFrame
    Gens_Init_State = CSV.File(events*"_gens$N.csv") |> DataFrame
    if ngi
        Gens_Init_State = natural_gas_interactions!(ps,Lines_Init_State,Gens_Init_State)
    end
    l_failures = Lines_Init_State.state;
    ps.branch.status[l_failures .== 0] .= 0;
    l_recovery_times = deepcopy(Lines_Init_State.recovery_time);
    # generator states and recovery times
    g_failures = ones(size(ps.gen,1));
    g_failures[1:length(Gens_Init_State.state)] = Gens_Init_State.state;
    ps.gen.status[g_failures .== 0] .= 0;
    ps.gen.Pg[g_failures .== 0] .= 0;
    g_recovery_times = zeros(size(ps.gen,1));
    g_recovery_times[1:length(Gens_Init_State.recovery_time)] = Gens_Init_State.recovery_time;
    ## run step 3
    dt = 60
    ti = 60*48;
    t0 = 10
    Restore = crisp_RLOPF_inter(ps,l_recovery_times,g_recovery_times,dt,
              ti,t0,gen_on,comm,nucp,ngi,crt)
              println(Restore)
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
        ## run step 4
        ResilienceTri = crisp_res(Restore);
        println(ResilienceTri)
    end
    ## make dataframes
    case_res = DataFrame(resilience = ResilienceTri);
    ## save data
    CSV.write("results/$out_folder", case_res);
end

function Resilience_interact(N,ps_folder,out_folder,events,dt,comm,nucp,ngi,crt;param_file = "")
    rng = MersenneTwister(100+N)
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
    # step 1
    Lines_Init_State = CSV.File(events*"_lines$N.csv") |> DataFrame
    Gens_Init_State = CSV.File(events*"_gens$N.csv") |> DataFrame
    if ngi
        Gens_Init_State = natural_gas_interactions!(ps,Lines_Init_State,Gens_Init_State)
    end
    l_failures = Lines_Init_State.state;
    ps.branch.status[l_failures .== 0] .= 0;
    l_recovery_times = Lines_Init_State.recovery_time;
    # generator states and recovery times
    g_failures = ones(size(ps.gen,1));
    g_failures[1:length(Gens_Init_State.state)] = Gens_Init_State.state;
    ps.gen.status[g_failures .== 0] .= 0;
    ps.gen.Pg[g_failures .== 0] .= 0;
    g_recovery_times = zeros(size(ps.gen,1));
    g_recovery_times[1:length(Gens_Init_State.recovery_time)] = Gens_Init_State.recovery_time;
    #=#check for islands
    subgraph = find_subgraphs(ps);
    M = Int64(findmax(subgraph)[1]);
    ps_islands = build_islands(subgraph,ps)
    for i in 1:M
        psi = ps_subset(ps,ps_islands[i]);
        crisp_dcpf_g1_s!(psi);
        # run lsopf
        #dt = 10; # first minute, affects ramp rate limits - rateB
        #crisp_lsopf_g1_s!(psi,dt);
        ps.gen.Pg[ps_islands[i].gen] = psi.gen.Pg
        #ps.storage.Ps[ps_islands[i].storage] = psi.storage.Ps
        #ps.storage.E[ps_islands[i].storage] = psi.storage.E
	    for sh in 1:length(psi.shunt.status)
		    if psi.shunt.status[sh] > 1
			     psi.shunt.status[sh] = 1.0
		    end
	    end
        ps.shunt.status[ps_islands[i].shunt] = psi.shunt.status
    end =#
    ## run step 3
    dt = 60
    ti = 60*48;
    t0 = 10
    Restore = crisp_Restoration_inter(ps,l_recovery_times,g_recovery_times,dt,
              ti,t0,gen_on,comm,nucp,ngi,crt)
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
        ## run step 4
        ResilienceTri = crisp_res(Restore);
        println(ResilienceTri)
    end
    ## make dataframes
    case_res = DataFrame(resilience = ResilienceTri);
    ## save data
    CSV.write("results/$out_folder", case_res);
end

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
    crisp_dcpf_g1_s!(ps)
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
            crisp_dcpf_g1_s!(psi);
            ps.gen.Pg[ps_islands[i].gen] = psi.gen.Pg
	    for sh in 1:length(psi.shunt.status)
		if psi.shunt.status[sh] > 1
			psi.shunt.status[sh] = 1.0
		end
	    end
            ps.shunt.status[ps_islands[i].shunt] = psi.shunt.status
        end
        ## run step 3
        dt = 60
        ti = 60*48;
        t0 = 10
        Restore = crisp_Restoration_var(ps,l_recovery_times,g_recovery_times,dt,ti,t0,gen_on)
        LoadShed0[iterat] = Restore.load_shed[2];
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

##CRISP_Rdist_gen_mh_varL+PV.jl
#Pulls from the distributions in the CRISP_initiate.jl file
#=
include("CRISP_initiate.jl")
include("CRISP_LSOPF_gen1.jl")
include("CRISP_RLOPF_mh_varLoad.jl")
include("CRISP_RT.jl")
include("CRISP_network_gen.jl")
=#
function Resilience_mh_var(Num,ps_folder,out_folder,dt;param_file = "")

    debug=1;
    tolerance1 = 10^(-6);
    ## Num = number of failure scenarios to run through
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
    #add columns to keep track of the time each generator is on or off
    if sum(names(ps.gen).==:time_on) == 0
        ps.gen.time_on = zeros(length(ps.gen.Pg));
    end
    if sum(names(ps.gen).==:time_off) == 0
        ps.gen.time_off = zeros(length(ps.gen.Pg));
        ps.gen.time_off[ps.gen.Pg.==0] .= ps.gen.minDownTimeHr[ps.gen.Pg.==0];
    end
    ps0 = deepcopy(ps);

    if isempty(param_file)
        # parameters of distributions for line outages and recovery times
        s_line = 2.56;#lines_dist[1];
        maxLinesOut = length(ps.branch.f); # => k in zipf distribution
        mu_line = 3.66;#lines_dist[2];
        sigma_line = 2.43;#lines_dist[3];
        lambda_gen = 1;
    end

    for iterat in 1:Num
        ps = deepcopy(ps0); # reset network to original state
        # step 1
        Lines_Init_State = line_state!(ps,s_line,maxLinesOut,mu_line,sigma_line)
        Gens_Init_State = gen_state!(ps,lambda_gen,mu_line,sigma_line)
        if debug==1
            outnow = (out_folder[1:end-4]);
            CSV.write("results"*outnow*"_lines.csv", Lines_Init_State)
            CSV.write("results"*outnow*"_gens.csv", Gens_Init_State)
        end
        l_failures = Lines_Init_State[:,1];
        ps.branch.status[l_failures .==0] .= 0;
        NumLinesOut = length(l_failures) - sum(l_failures);
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
            # run lsopf
            dt = 10;
            crisp_lsopf_g1_s!(psi,dt);
            ps.gen.Pg[ps_islands[i].gen] = psi.gen.Pg
            ps.storage.Ps[ps_islands[i].storage] = psi.storage.Ps
            #ps.storage.E[ps_islands[i].storage] = psi.storage.E
            ps.shunt.status[ps_islands[i].shunt] = psi.shunt.status
        end
        println(iterat)
        @assert total>=sum(ps.shunt.P .* ps.shunt.status)
        @assert 2M*tolerance1 >= abs(sum(ps.shunt.P .*ps.shunt.status)-sum(ps.gen.Pg)-sum(ps.storage.Ps))
        tolerance = 10^(-10);
        LoadShed0[iterat] = total-sum(ps.shunt.P .* ps.shunt.status);
        ## run step 3
        dt = 15
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


##CRISP_Rdist_gen_movh.jl
#=
include("CRISP_initiate.jl")
include("CRISP_LSOPF_gen1.jl")
include("CRISP_RLOPF_movh.jl")
include("CRISP_RT.jl")
include("CRISP_network_gen.jl")
=#
function Res_dist_mh(Num,ps_folder,out_folder,dt;param_file = "")
    debug=1;
    tolerance1 = 10^(-6);
    ## Num = number of failure scenarios to run through
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
    total = sum(ps.shunt[:P]);
    Pd_max = deepcopy(ps.shunt[:P]);
    ps0 = deepcopy(ps);

    if isempty(param_file)
        # parameters of distributions for line outages and recovery times
        s_line = 2.56;#lines_dist[1];
        maxLinesOut = length(ps.branch.f); # => k in zipf distribution
        mu_line = 3.66;#lines_dist[2];
        sigma_line = 2.43;#lines_dist[3];
        lambda_gen = 1;
        #pick start up times for gens:
        gen_startup = 90 .*ones(length(ps.gen.bus));
    end

    for iterat in 1:Num
        ps = deepcopy(ps0); # reset network to original state
        # step 1
        Lines_Init_State = line_state!(ps,s_line,maxLinesOut,mu_line,sigma_line)
        Gens_Init_State = gen_state!(ps,lambda_gen,mu_line,sigma_line)
        if debug==1
            outnow = (out_folder[1:end-4]);
            CSV.write("results$outnow IC$iterat lines.csv", Lines_Init_State)
            CSV.write("results$outnow IC$iterat gens.csv", Gens_Init_State)
        end
        l_failures = Lines_Init_State[:,1];
        ps.branch.status[l_failures .==0] .= 0;
        NumLinesOut = length(l_failures) - sum(l_failures);
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
            # run lsopf
            dt = 10;
            crisp_lsopf_g1_s!(psi,dt);
            ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
            ps.storage[ps_islands[i].storage,:Ps] = psi.storage.Ps
            ps.storage[ps_islands[i].storage,:E] = psi.storage.E
            ps.shunt[ps_islands[i].shunt,:status] = psi.shunt.status
        end
        println(iterat)
        @assert total>=sum(ps.shunt.P .* ps.shunt.status)
        @assert 2M*tolerance1 >= abs(sum(ps.shunt.P .*ps.shunt.status)-sum(ps.gen.Pg)-sum(ps.storage.Ps))
        tolerance = 10^(-10);
        LoadShed0[iterat] = total-sum(ps.shunt.P .* ps.shunt.status);
        ## run step 3
        dt = 1
        ti = dt;#10
        t0 = 10
        #crisp_mh_rlopf!(ps,dt,time)
        Restore = crisp_Restore_mh(ps,l_recovery_times,g_recovery_times,dt,ti,t0)
        if debug==1
            outnow = (out_folder[1:end-4]);
            CSV.write("results$outnow IC$iterat restore.csv", Restore)
        end
        #Restore = RLSOPF_g_s!(ps,dt,state,gens_state,recovery_times,gens_recovery_time,Pd_max)# data frame [times, load shed in cost per hour]
        ###find the time to restore the grid to 99.9% load served
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
        end
    end
    case_res = DataFrame(resilience = ResilienceTri[:,1]);
    #case_lines = DataFrame(lines_out = NumLinesOut[:,1], no_load_shed_time = LoadServedTime[:,1], initial_load_shed = LoadShed0[:,1]);#, full_rest_time = MaxRestorationTime[:,1]);

    ## save data
    CSV.write("results\\$out_folder", case_res);
    #CSV.write("results\\$(out_folder[1:end-4])_lines.csv", case_lines);
end

##CRISP_Rdist_gen_stor1.jl
#=
include("CRISP_initiate.jl")
include("CRISP_LSOPF_gen_stor.jl")
include("CRISP_RLOPF_gen_stor.jl")
include("CRISP_RT.jl")
include("CRISP_network_gen.jl")
=#
function Res_dist_gen_stor1(Num,ps_folder,out_folder,dt;param_file = "")
    debug=0;
    ## Num = number of failure scenarios to run through
    # initialize vector of costs from events
    NumLinesOut = Array{Float64}(undef,Num);
    LoadShed0 =  Array{Float64}(undef,Num);
    MaxRestorationTime = Array{Float64}(undef,Num);
    LoadServedTime = Array{Float64}(undef,Num);
    ResilienceTri = Array{Float64}(undef,Num);
    ## load the case data
    ps = import_ps("$ps_folder")
    ps.shunt = ps.shunt[ps.shunt.P .!=0,:];
    #add columns to keep track of the time each generator is on or off
    if sum(names(ps.gen).==:time_on) == 0
        ps.gen.time_on = zeros(length(ps.gen.Pg));
    end
    if sum(names(ps.gen).==:time_off) == 0
        ps.gen.time_off = zeros(length(ps.gen.Pg));
    end
    crisp_opf_initiate!(ps,0.1)
    total = sum(ps.shunt[:P]);
    Pd_max = deepcopy(ps.shunt[:P]);
    ps0 = deepcopy(ps);
    if isempty(param_file)
        # parameters of distributions for line outages and recovery times
        s_line = 2.56;#lines_dist[1];
        maxLinesOut = length(ps.branch.f); # => k in zipf distribution
        mu_line = 3.66;#lines_dist[2];
        sigma_line = 2.43;#lines_dist[3];
        lambda_gen = 1;
        #pick start up times for gens:
        gen_startup = 90 .*ones(length(ps.gen.bus));
    end

    for iterat in 1:Num
        ps = deepcopy(ps0); # reset network to original state
        # step 1
        Lines_Init_State = line_state!(ps,s_line,maxLinesOut,mu_line,sigma_line)
        Gens_Init_State = gen_state!(ps,lambda_gen,mu_line,sigma_line)
        if debug==1
            outnow = (out_folder[1:end-4]);
            CSV.write("results$outnow IC$iterat lines.csv", Lines_Init_State)
            CSV.write("results$outnow IC$iterat gens.csv", Gens_Init_State)
        end
        l_failures = Lines_Init_State[:,1];
        ps.branch.status[l_failures .==0] .= 0;
        NumLinesOut = length(l_failures) - sum(l_failures);
        l_recovery_times = Lines_Init_State[:,2];
        # generator states and recovery times
        g_failures = Gens_Init_State[:,1];
        ps.gen.status[g_failures.==0] .= 0;
        ps.gen.Pg[g_failures.==0] .= 0;
        g_recovery_times = Gens_Init_State[:,2];
        #check for islands
        subgraph = find_subgraphs(ps);
        M = Int64(findmax(subgraph)[1]);
        ps_islands = build_islands(subgraph,ps);
        for i in 1:M
            psi = ps_subset(ps,ps_islands[i]);
            # run lsopf
            crisp_lsopf_g_s!(psi,dt);
            ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
            ps.shunt[ps_islands[i].shunt,:status] = psi.shunt.status
            ps.storage[ps_islands[i].storage,:E] = psi.storage.E
            ps.storage[ps_islands[i].storage,:Ps] = psi.storage.Ps
        end
        @assert 10^(-6)>=abs(sum(ps.shunt.P .*ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg))
        @assert total>=sum(ps.shunt.P)
        println(iterat)
        @assert total>=sum(ps.shunt.P .*ps.shunt.status)
        @assert 10^(-6)>=abs(sum(ps.shunt.P .*ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg))
        LoadShed0[iterat] = total-sum(ps.shunt.P.*ps.shunt.status);
        ## run step 3
        dt = 1
        ti = dt;#10
        t0 = 10
        #crisp_mh_rlopf!(ps,dt,time)
        Restore = crisp_Restore(ps,l_recovery_times,g_recovery_times,dt,ti,t0)
        if debug==1
            CSV.write("results$outnow IC$iterat restore.csv", Restore)
        end
        ###find the time to restore the grid to 99.9% load served
        K = abs.(Restore.perc_load_served .- 1) .<= 0.001;
        K[1] = false;
        if isempty( Restore.time[K]) && (size(Restore)[1] > 1)
            error("Never solved to full restoration.")
        elseif isempty( Restore.time[K]) && (size(Restore)[1] == 1)
            LoadServedTime[iterat] = 0.0;
        else
            R = Restore.time[K];
            LoadServedTime[iterat] = R[1] - Restore.time[1];
            ## run step 4
            ResilienceTri[iterat] = crisp_res(Restore);
        end
    end
    case_res = DataFrame(resilience = ResilienceTri[:,1]);
    case_lines = DataFrame(lines_out = NumLinesOut, no_load_shed_time = LoadServedTime, initial_load_shed = LoadShed0, full_rest_time = MaxRestorationTime);
    ## save data
    CSV.write("results\\$out_folder", case_res);
    CSV.write("results\\$(out_folder[1:end-4])_info.csv", case_lines);
end


##CRISP_gen_stor.jl
#=
include("CRISP_initiate.jl")
include("CRISP_LSOPF_gen.jl")
include("CRISP_RLSOPF_gen.jl")
include("CRISP_RT.jl")
include("CRISP_network_gen.jl")
=#
function Res_dist_gen_stor(Num,ps_folder,out_folder,dt;param_file = "")
    debug=1;
    ## Num = number of failure scenarios to run through
    # initialize vector of costs from events
    NumLinesOut = Array{Float64}(undef,Num,1);
    LoadShed0 =  Array{Float64}(undef,Num,1);
    MaxRestorationTime = Array{Float64}(undef,Num,1);
    LoadServedTime = Array{Float64}(undef,Num,1);
    ResilienceTri = Array{Float64}(undef,Num,1);
    ## load the case data
    ps = import_ps("$ps_folder")
    crisp_dcpf_g!(ps)
    total = sum(ps.shunt[:P]);
    Pd_max = deepcopy(ps.shunt[:P]);
    ps0 = deepcopy(ps);

    if isempty(param_file)
        # parameters of distributions for line outages and recovery times
        s_line = 2.56;#lines_dist[1];
        maxLinesOut = length(ps.branch.f); # => k in zipf distribution
        mu_line = 3.66;#lines_dist[2];
        sigma_line = 2.43;#lines_dist[3];
        lambda_gen = 1;
        #pick start up times for gens:
        gen_startup = 90 .*ones(length(ps.gen.bus));
    end

    for iterat in 1:Num
        ps = deepcopy(ps0); # reset network to original state
        # step 1
        Lines_Init_State = line_state!(ps,s_line,maxLinesOut,mu_line,sigma_line)
        Gens_Init_State = gen_state!(ps,lambda_gen,mu_line,sigma_line)
        if debug==1
            outnow = (out_folder[1:end-4]);
            CSV.write("results$outnow IC$iterat lines.csv", Lines_Init_State)
            CSV.write("results$outnow IC$iterat gens.csv", Gens_Init_State)
        end
        state = Lines_Init_State[:,1];
        NumLinesOut[iterat] = length(state) - sum(state)
        recovery_times = Lines_Init_State[:,2];
        MaxRestorationTime[iterat] = maximum(recovery_times);
        failures = state;
        # generator states and recovery times
        gens_state = Gens_Init_State[:,1];
        gens_recovery_time = Gens_Init_State[:,2];
        #check for islands
        subgraph = find_subgraphs(ps);
        M = Int64(findmax(subgraph)[1]);
        if M>1
            ps_islands = build_islands(subgraph,ps);
            for i in 1:M
                psi = ps_subset(ps,ps_islands[i]);
                # run the dcpf
                crisp_dcpf_g_s!(psi);
                # run lsopf
                crisp_lsopf_g_s!(psi,dt);
                ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
                ps.shunt[ps_islands[i].shunt,:P] = psi.shunt.P
                ps.storage[ps_islands[i].storage,:E] = psi.storage.E
                ps.storage[ps_islands[i].storage,:Ps] = psi.storage.Ps
                crisp_dcpf_g_s!(psi);
            end
            @assert 10^(-6)>=abs(sum(ps.shunt.P)-sum(ps.gen.Pg))
            @assert total>=sum(ps.shunt.P)
        else
            crisp_dcpf_g_s!(ps);
            crisp_lsopf_g_s!(ps,dt);
            crisp_dcpf_g_s!(ps);
        end
        println(iterat)
        @assert total>=sum(ps.shunt.P)
        @assert 10^(-6)>=abs(sum(ps.shunt.P)-sum(ps.gen.Pg))
        tolerance = 10^(-10);
        if (abs(total-sum(ps.shunt.P)) <= tolerance)
            ResilienceTri[iterat] = 0;
            LoadServedTime[iterat] = 0;
            LoadShed0[iterat] = 0;
        else
            LoadShed0[iterat] = total-sum(ps.shunt.P);
            ## run step 3
            Restore = RLSOPF_g_s!(ps,dt,state,gens_state,recovery_times,gens_recovery_time,Pd_max)# data frame [times, load shed in cost per hour]
            ###find the time to restore the grid to 99.9% load served
            K = abs.(Restore.perc_load_served .- 1) .<= 0.001;
            K[1] = false;
            if isempty( Restore.time[K])
                error("Never solved to full restoration.")
            else
                R = Restore.time[K];
                LoadServedTime[iterat] = R[1] - Restore.time[2];
                ## run step 4
                ResilienceTri[iterat] = crisp_res(Restore);
            end
        end
    end
    case_res = DataFrame(resilience = ResilienceTri[:,1]);
    case_lines = DataFrame(lines_out = NumLinesOut[:,1], no_load_shed_time = LoadServedTime[:,1], initial_load_shed = LoadShed0[:,1], full_rest_time = MaxRestorationTime[:,1]);

    ## save data
    CSV.write("results\\$out_folder", case_res);
    CSV.write("results\\$(out_folder[1:end-4])_lines.csv", case_lines);
end

## CRISP_Rdist_gen.jl
#=
include("CRISP_initiate.jl")
include("CRISP_LSOPF_gen.jl")
include("CRISP_RLSOPF_gen.jl")
include("CRISP_RT.jl")
include("CRISP_network_gen.jl")
=#
function Res_dist_gen(Num,ps_folder,out_folder;param_file = "")
    debug=0;
    ## Num = number of failure scenarios to run through
    # initialize vector of costs from events
    NumLinesOut = Array{Float64}(undef,Num,1);
    LoadShed0 =  Array{Float64}(undef,Num,1);
    MaxRestorationTime = Array{Float64}(undef,Num,1);
    LoadServedTime = Array{Float64}(undef,Num,1);
    ResilienceTri = Array{Float64}(undef,Num,1);
    ## load the case data
    ps = import_ps("$ps_folder")
    crisp_dcpf_g!(ps)
    total = sum(ps.shunt[:P]);
    Pd_max = deepcopy(ps.shunt[:P]);
    ps0 = deepcopy(ps);

    if isempty(param_file)
        # parameters of distributions for line outages and recovery times
        s_line = 2.56;#lines_dist[1];
        maxLinesOut = length(ps.branch.f); # => k in zipf distribution
        mu_line = 3.66;#lines_dist[2];
        sigma_line = 2.43;#lines_dist[3];
        lambda_gen = 1;
        #pick start up times for gens:
        gen_startup = 90 .*ones(length(ps.gen.bus));
    end

    for iterat in 1:Num
        ps = deepcopy(ps0); # reset network to original state
        # step 1
        Lines_Init_State = line_state!(ps,s_line,maxLinesOut,mu_line,sigma_line)
        Gens_Init_State = gen_state!(ps,lambda_gen,mu_line,sigma_line)
        if debug==1
            outnow = (out_folder[1:end-4]);
            CSV.write("results$outnow IC$iterat lines.csv", Lines_Init_State)
            CSV.write("results$outnow IC$iterat gens.csv", Gens_Init_State)
        end
        state = Lines_Init_State[:,1];
        NumLinesOut[iterat] = length(state) - sum(state)
        recovery_times = Lines_Init_State[:,2];
        MaxRestorationTime[iterat] = maximum(recovery_times);
        failures = state;
        # generator states and recovery times
        gens_state = Gens_Init_State[:,1];
        gens_recovery_time = Gens_Init_State[:,2];
        #check for islands
        subgraph = find_subgraphs(ps);
        M = Int64(findmax(subgraph)[1]);
        if M>1
            ps_islands = build_islands(subgraph,ps);
            for i in 1:M
                psi = ps_subset(ps,ps_islands[i]);
                # run the dcpf
                crisp_dcpf_g!(psi);
                # run lsopf
                crisp_lsopf_g!(psi);
                ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
                ps.shunt[ps_islands[i].shunt,:P] = psi.shunt.P
                crisp_dcpf_g!(psi);
            end
            @assert 10^(-6)>=abs(sum(ps.shunt.P)-sum(ps.gen.Pg))
            @assert total>=sum(ps.shunt.P)
        else
            crisp_dcpf_g!(ps);
            crisp_lsopf_g!(ps);
            crisp_dcpf_g!(ps);
        end
        println(iterat)
        @assert total>=sum(ps.shunt.P)
        @assert 10^(-6)>=abs(sum(ps.shunt.P)-sum(ps.gen.Pg))
        tolerance = 10^(-10);
        if (abs(total-sum(ps.shunt.P)) <= tolerance)
            ResilienceTri[iterat] = 0;
            LoadServedTime[iterat] = 0;
            LoadShed0[iterat] = 0;
        else
            LoadShed0[iterat] = total-sum(ps.shunt.P);
            ## run step 3
            Restore = RLSOPF_g!(ps,state,gens_state,recovery_times,gens_recovery_time,Pd_max)# data frame [times, load shed in cost per hour]
            ###find the time to restore the grid to 99.9% load served
            K = abs.(Restore.perc_load_served .- 1) .<= 0.001;
            K[1] = false;
            if isempty( Restore.time[K])
                error("Never solved to full restoration.")
            else
                R = Restore.time[K];
                LoadServedTime[iterat] = R[1] - Restore.time[2];
                ## run step 4
                ResilienceTri[iterat] = crisp_res(Restore);
            end
        end
    end
    case_res = DataFrame(resilience = ResilienceTri[:,1]);
    case_lines = DataFrame(lines_out = NumLinesOut[:,1], no_load_shed_time = LoadServedTime[:,1], initial_load_shed = LoadShed0[:,1], full_rest_time = MaxRestorationTime[:,1]);

    ## save data
    CSV.write("results\\$out_folder", case_res);
    CSV.write("results\\$(out_folder[1:end-4])_lines.csv", case_lines);
end

##CRISP_Rdist_gen_orig.jl
#=
include("CRISP_initiate.jl")
include("CRISP_LSOPF_gen.jl")
include("CRISP_RLSOPF_gen_orig.jl")
include("CRISP_RT.jl")
include("CRISP_network.jl")
=#
function Res_dist_gen_orig(Num,ps_folder,out_folder;param_file = "")
    debug=0;
    ## Num = number of failure scenarios to run through
    # initialize vector of costs from events
    NumLinesOut = Array{Float64}(undef,Num,1);
    LoadShed0 =  Array{Float64}(undef,Num,1);
    MaxRestorationTime = Array{Float64}(undef,Num,1);
    LoadServedTime = Array{Float64}(undef,Num,1);
    ResilienceTri = Array{Float64}(undef,Num,1);
    ## load the case data
    ps = import_ps("$ps_folder")
    crisp_dcpf!(ps)
    total = sum(ps.shunt[:P]);
    Pd_max = deepcopy(ps.shunt[:P]);
    ps0 = deepcopy(ps);

    if isempty(param_file)
        # parameters of distributions for line outages and recovery times
        s_line = 2.56;#lines_dist[1];
        maxLinesOut = length(ps.branch.f); # => k in zipf distribution
        mu_line = 3.66;#lines_dist[2];
        sigma_line = 2.43;#lines_dist[3];
        lambda_gen = 1;
        #pick start up times for gens:
        gen_startup = 90 .*ones(length(ps.gen.bus));
    end

    for iterat in 1:Num
        ps = deepcopy(ps0); # reset network to original state
        # step 1
        Lines_Init_State = line_state!(ps,s_line,maxLinesOut,mu_line,sigma_line)
        Gens_Init_State = gen_state!(ps,lambda_gen,mu_line,sigma_line)
        if debug==1
            outnow = (out_folder[1:end-4]);
            CSV.write("results$outnow IC$iterat.csv", Lines_Init_State)
        end
        state = Lines_Init_State[:,1];
        NumLinesOut[iterat] = length(state) - sum(state)
        recovery_times = Lines_Init_State[:,2];
        MaxRestorationTime[iterat] = maximum(recovery_times);
        failures = state;
        # generator states and recovery times
        gens_state = Gens_Init_State[:,1];
        gens_recovery_time = Gens_Init_State[:,2];
        #check for islands
        subgraph = find_subgraphs(ps);
        M = Int64(findmax(subgraph)[1]);
        if M>1
            ps_islands = build_islands(subgraph,ps);
            for i in 1:M
                psi = ps_subset(ps,ps_islands[i]);
                # run the dcpf
                crisp_dcpf_g!(psi);
                # run lsopf
                crisp_lsopf_g!(psi);
                ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
                ps.shunt[ps_islands[i].shunt,:P] = psi.shunt.P
                crisp_dcpf_g!(psi);
            end
            @assert 10^(-6)>=abs(sum(ps.shunt.P)-sum(ps.gen.Pg))
            @assert total>=sum(ps.shunt.P)
        else
            crisp_dcpf_g!(ps);
            crisp_lsopf_g!(ps);
            crisp_dcpf_g!(ps);
        end
        println(iterat)
        @assert total>=sum(ps.shunt.P)
        @assert 10^(-6)>=abs(sum(ps.shunt.P)-sum(ps.gen.Pg))
        tolerance = 10^(-10);
        if (abs(total-sum(ps.shunt.P)) <= tolerance)
            ResilienceTri[iterat] = 0;
            LoadServedTime[iterat] = 0;
            LoadShed0[iterat] = 0;
        else
            LoadShed0[iterat] = total-sum(ps.shunt.P);
            ## run step 3
            Restore = RLSOPF_g_orig!(ps,state,gens_state,recovery_times,gens_recovery_time,gen_startup,Pd_max)# data frame [times, load shed in cost per hour]
            K = abs.(Restore.perc_load_served .- 1) .<= 0.001;
            K[1] = false;
            R = Restore.time[K];
            LoadServedTime[iterat] = R[1] - Restore.time[2];
            ## run step 4
            ResilienceTri[iterat] = crisp_res(Restore);
        end
    end
    case_res = DataFrame(resilience = ResilienceTri[:,1]);
    case_lines = DataFrame(lines_out = NumLinesOut[:,1], no_load_shed_time = LoadServedTime[:,1], initial_load_shed = LoadShed0[:,1], full_rest_time = MaxRestorationTime[:,1]);

    ## save data
    CSV.write("results\\$out_folder", case_res);
    CSV.write("results\\$(out_folder[1:end-4])_lines.csv", case_lines);
end

## Original CRISP_Rdist.jl
function Res_dist(Num,ps_folder,out_folder;param_file = "")
    debug=1;
    ## Num = number of failure scenarios to run through
    # initialize vector of costs from events
    NumLinesOut = Array{Float64}(undef,Num,1);
    LoadShed0 =  Array{Float64}(undef,Num,1);
    MaxRestorationTime = Array{Float64}(undef,Num,1);
    LoadServedTime = Array{Float64}(undef,Num,1);
    ResilienceTri = Array{Float64}(undef,Num,1);
    ## load the case data
    ps = import_ps0("$ps_folder")
    ps.shunt.P .= 2*ps.shunt.P;
    crisp_dcpf!(ps)
    total = sum(ps.shunt[:P]);
    Pd_max = deepcopy(ps.shunt[:P]);
    ps0 = deepcopy(ps);

    if isempty(param_file)
        # parameters of distributions for line outages and recovery times
        s_line = 2.56;#lines_dist[1];
        maxLinesOut = length(ps.branch.f); # => k in zipf distribution
        mu_line = 3.66;#lines_dist[2];
        sigma_line = 2.43;#lines_dist[3];
    end

    for iterat in 1:Num
        ps = deepcopy(ps0); # reset network to original state
        # step 1
        Lines_Init_State = line_state!(ps,s_line,maxLinesOut,mu_line,sigma_line)
        if debug==1
            outnow = (out_folder[1:end-4]);
            CSV.write("results$outnow IC$iterat.csv", Lines_Init_State)
        end
        state = Lines_Init_State[:,1];
        NumLinesOut[iterat] = length(state) - sum(state)
        recovery_times = Lines_Init_State[:,2];
        MaxRestorationTime[iterat] = maximum(recovery_times);
        failures = state;
        if sum(failures.==0)==0
            ResilienceTri[iterat] = 0;
            LoadServedTime[iterat] = 0;
            LoadShed0[iterat] = 0;
            println(iterat)
        else
            #check for islands
            subgraph = find_subgraphs0(ps);
            M = Int64(findmax(subgraph)[1]);
            if M>1
                ps_islands = build_islands0(subgraph,ps);
                for i in 1:M
                    psi = ps_subset0(ps,ps_islands[i]);
                    # run the dcpf
                    crisp_dcpf!(psi);
                    # run lsopf
                    crisp_lsopf1!(psi);
                    ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
                    ps.shunt[ps_islands[i].shunt,:P] = psi.shunt.P
                    crisp_dcpf!(psi);
                end
                    @assert 10^(-6)>=abs(sum(ps.shunt.P)-sum(ps.gen.Pg))
                    @assert total>=sum(ps.shunt.P)
            else
                crisp_dcpf!(ps);
                crisp_lsopf1!(ps);
                crisp_dcpf!(ps);
            end
            println(iterat)
            @assert total>=sum(ps.shunt.P)
            @assert 10^(-6)>=abs(sum(ps.shunt.P)-sum(ps.gen.Pg))
            tolerance = 10^(-10);
            if (abs(total-sum(ps.shunt.P)) <= tolerance)
                ResilienceTri[iterat] = 0;
                LoadServedTime[iterat] = 0;
                LoadShed0[iterat] = 0;
            else
                LoadShed0[iterat] = total-sum(ps.shunt.P);
                ## run step 3
                Restore = RLSOPF!(total,ps,failures,recovery_times,Pd_max);#,load_cost) # data frame [times, load shed in cost per hour]
                println(Restore)
                K = abs.(Restore.perc_load_served .- 1) .<= 0.001;
                K[1] = false;
                R = Restore.time[K];
                LoadServedTime[iterat] = R[1] - Restore.time[2];
                ## run step 4
                ResilienceTri[iterat] = crisp_res0(Restore);
            end
        end
    end
    case_res = DataFrame(resilience = ResilienceTri[:,1]);
    case_lines = DataFrame(lines_out = NumLinesOut[:,1], no_load_shed_time = LoadServedTime[:,1], initial_load_shed = LoadShed0[:,1], full_rest_time = MaxRestorationTime[:,1]);

    ## save data
    CSV.write("results\\$out_folder", case_res);
    CSV.write("results\\$(out_folder[1:end-4])_lines.csv", case_lines);
end

function Res_dist_test2(Num,ps_folder,out_folder;param_file = "")
    debug=0;
    ## Num = number of failure scenarios to run through
    # initialize vector of costs from events
    ResilienceTri = Array{Float64}(undef,Num,1);
    ## load the case data
    ps = import_ps0("$ps_folder")
    crisp_dcpf!(ps)
    total = sum(ps.shunt[:P]);
    Pd_max = deepcopy(ps.shunt[:P]);
    ps0 = deepcopy(ps);

    if isempty(param_file)
        # parameters of distributions for line outages and recovery times
        s_line = 2.56;#lines_dist[1];
        maxLinesOut = length(ps.branch.f); # => k in zipf distribution
        mu_line = 3.66;#lines_dist[2];
        sigma_line = 2.43;#lines_dist[3];
    end

    for iterat in 1:Num
        ps = deepcopy(ps0); # reset network to original state
        # step 1
        Lines_Init_State = line_state!(ps,s_line,maxLinesOut,mu_line,sigma_line)
        if debug==1
            CSV.write("results\\$out_folder$iterat.csv", Lines_Init_State)
        end
        state = Lines_Init_State[:,1];
        recovery_times = Lines_Init_State[:,2];
        failures = state;
        if sum(failures.==0)==0
            ResilienceTri[iterat] = 0;
            println(iterat)
        else
            #check for islands
            subgraph = find_subgraphs0(ps);
            M = Int64(findmax(subgraph)[1]);
            ps_islands = build_islands0(subgraph,ps);
            for i in 1:M
                psi = ps_subset0(ps,ps_islands[i]);
                # run the dcpf
                crisp_dcpf!(psi);
                # run lsopf
                crisp_lsopf!(psi);
                ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
                ps.shunt[ps_islands[i].shunt,:P] = psi.shunt.P
                crisp_dcpf!(psi);
            end
            @assert 10^(-6)>=abs(sum(ps.shunt.P)-sum(ps.gen.Pg))
            @assert total>=sum(ps.shunt.P)
            println(iterat)
            @assert total>=sum(ps.shunt.P)
            @assert 10^(-6)>=abs(sum(ps.shunt.P)-sum(ps.gen.Pg))
            tolerance = 10^(-10);
            if (abs(total-sum(ps.shunt.P)) <= tolerance)
                ResilienceTri[iterat] = 0;
            else
                ## run step 3
                Restore = RLSOPF!(total,ps,failures,recovery_times,Pd_max);#,load_cost) # data frame [times, load shed in cost per hour]
                ## run step 4
                ResilienceTri[iterat] = crisp_res0(Restore);
            end
        end
    end
    case_res = DataFrame(resilience = ResilienceTri[:,1]);
    ## save data
    CSV.write("results\\$out_folder", case_res);
end

## CRISP_Rdist_0Out.jl

function Res_dist_0O(Num,ps_folder,out_folder;param_file = "")
    debug=0;
    ## Num = number of failure scenarios to run through
    # initialize vector of costs from events
    ResilienceTri = Array{Float64}(undef,Num,1);
    ## load the case data
    ps = import_ps0("$ps_folder")
    crisp_dcpf!(ps)
    total = sum(ps.shunt[:P]);
    Pd_max = deepcopy(ps.shunt[:P]);
    ps0 = deepcopy(ps);

    if isempty(param_file)
        # parameters of distributions for line outages and recovery times
        s_line = 2.56;#lines_dist[1];
        maxLinesOut = length(ps.branch.f); # => k in zipf distribution
        mu_line = 3.66;#lines_dist[2];
        sigma_line = 2.43;#lines_dist[3];
    end

    for iterat in 1:Num
        ps = deepcopy(ps0); # reset network to original state
        # step 1
        Lines_Init_State = line_state2!(ps,s_line,maxLinesOut,mu_line,sigma_line)
        if debug==1
            CSV.write("results\\$out_folder$iterat.csv", Lines_Init_State)
        end
        state = Lines_Init_State[:,1];
        recovery_times = Lines_Init_State[:,2];
        failures = state;
        if sum(failures.==0)==0
            ResilienceTri[iterat] = 0;
            println(iterat)
        else
            #check for islands
            subgraph = find_subgraphs0(ps);
            M = Int64(findmax(subgraph)[1]);
            if M>1
                ps_islands = build_islands0(subgraph,ps);
                for i in 1:M
                    psi = ps_subset0(ps,ps_islands[i]);
                    # run the dcpf
                    crisp_dcpf!(psi);
                    # run lsopf
                    crisp_lsopf!(psi);
                    ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
                    ps.shunt[ps_islands[i].shunt,:P] = psi.shunt.P
                    crisp_dcpf!(psi);
                end
                    @assert 10^(-6)>=abs(sum(ps.shunt.P)-sum(ps.gen.Pg))
                    @assert total>=sum(ps.shunt.P)
            else
                crisp_dcpf!(ps);
                crisp_lsopf!(ps);
                crisp_dcpf!(ps);
            end
            println(iterat)
            @assert total>=sum(ps.shunt.P)
            @assert 10^(-6)>=abs(sum(ps.shunt.P)-sum(ps.gen.Pg))
            tolerance = 10^(-10);
            if (abs(total-sum(ps.shunt.P)) <= tolerance)
                ResilienceTri[iterat] = 0;
            else
                ## run step 3
                Restore = RLSOPF!(total,ps,failures,recovery_times,Pd_max);#,load_cost) # data frame [times, load shed in cost per hour]
                ## run step 4
                ResilienceTri[iterat] = crisp_res0(Restore);
            end
        end
    end
    case_res = DataFrame(resilience = ResilienceTri[:,1]);
    ## save data
    CSV.write("results\\$out_folder", case_res);

end
