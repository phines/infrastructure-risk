using CSV
include("CRISP_initiate.jl")
include("CRISP_LSOPF_gen.jl")
include("CRISP_RLSOPF_gen.jl")
include("CRISP_RT.jl")
include("CRISP_network_gen.jl")

function Res_dist_gen_stor(Num,ps_folder,out_folder,dt;param_file = "")
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


#=
#for debugging
include("src\\CRISP_initiate.jl")
include("src\\CRISP_LSOPF_gen.jl")
include("src\\CRISP_RLSOPF_gen.jl")
include("src\\CRISP_RT.jl")
include("src\\CRISP_network_gen.jl")
=#