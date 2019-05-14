using CSV
include("CRISP_initiate.jl")
include("CRISP_LSOPF.jl")
include("CRISP_RLSOPF.jl")
include("CRISP_RT.jl")
include("CRISP_network.jl")

function Res_dist(Num,ps_folder,out_folder;param_file = "")
    debug=1;
    ## Num = number of failure scenarios to run through
    # initialize vector of costs from events
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
        recovery_times = Lines_Init_State[:,2];
        failures = state;
        if sum(failures.==0)==0
            ResilienceTri[iterat] = 0;
            println(iterat)
        else
            #check for islands
            subgraph = find_subgraphs(ps);
            M = Int64(findmax(subgraph)[1]);
            if M>1
                ps_islands = build_islands(subgraph,ps);
                for i in 1:M
                    psi = ps_subset(ps,ps_islands[i]);
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
                ResilienceTri[iterat] = crisp_res(Restore);
            end
        end
    end
    case_res = DataFrame(resilience = ResilienceTri[:,1]);
    ## save data
    CSV.write("results\\$out_folder", case_res);

end

function Res_dist_test2(Num,ps_folder,out_folder;param_file = "")
    debug=0;
    ## Num = number of failure scenarios to run through
    # initialize vector of costs from events
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
            subgraph = find_subgraphs(ps);
            M = Int64(findmax(subgraph)[1]);
            #if M>1
                ps_islands = build_islands(subgraph,ps);
                for i in 1:M
                    psi = ps_subset(ps,ps_islands[i]);
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
            #else
                #crisp_dcpf!(ps);
                #crisp_lsopf!(ps);
                #crisp_dcpf!(ps);
            #end
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
                ResilienceTri[iterat] = crisp_res(Restore);
            end
        end
    end
    case_res = DataFrame(resilience = ResilienceTri[:,1]);
    ## save data
    CSV.write("results\\$out_folder", case_res);

end

#for debugging
#include("src\\CRISP_initiate.jl")
#include("src\\CRISP_LSOPF_tests.jl")
#include("src\\CRISP_RLSOPF_test.jl")
#include("src\\CRISP_RT.jl")
#include("src\\CRISP_network.jl")
