using CSV
include("CRISP_initiate.jl")
include("CRISP_LSOPF_gen_stor.jl")
include("CRISP_RLOPF_gen_stor.jl")
include("CRISP_RT.jl")
include("CRISP_network_gen.jl")

function Res_dist_gen_stor(Num,ps_folder,out_folder,dt;param_file = "")
    debug=1;
    outnow = (out_folder[1:end-4])
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
            CSV.write("results"*outnow*"_restore.csv", Restore)
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
    CSV.write("results"*out_folder, case_res);
    CSV.write("results"*(out_folder[1:end-4])*"_info.csv", case_lines);
end


#=
#for debugging
include("src\\CRISP_initiate.jl")
include("src\\CRISP_LSOPF_gen_stor.jl")
include("src\\CRISP_RLSOPF_gen_stor.jl")
include("src\\CRISP_RT.jl")
include("src\\CRISP_network_gen.jl")
=#
