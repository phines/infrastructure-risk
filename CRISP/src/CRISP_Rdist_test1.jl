using CSV
include("CRISP_initiate.jl")
include("CRISP_LSOPF_tests.jl")
include("CRISP_RLSOPF_test.jl")
include("CRISP_RT.jl")
include("CRISP_network.jl")

function Res_dist_test1(Num,ps_folder,out_folder;param_file = "")
    ## Num = number of failure scenarios to run through
    # initialize vector of costs from events
    ResilienceTri = Array{Float64}(undef,Num,1);
    ## load the case data
    ps = import_ps("$ps_folder")
    crisp_dcpf1!(ps)
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
        state = Lines_Init_State[:,1];
        recovery_times = Lines_Init_State[:,2];
        failures = state;

        #check for islands
        subgraph = find_subgraphs(ps);
        M = Int64(findmax(subgraph)[1]);
        ps_islands = build_islands(subgraph,ps);
        for i in 1:M
            psi = ps_subset(ps,ps_islands[i]);
            # run the dcpf
            crisp_dcpf1!(psi);
            # run lsopf
            crisp_lsopf1!(psi);
            ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
            ps.shunt[ps_islands[i].shunt,:P] = psi.shunt.P
            crisp_dcpf1!(psi);
        end

        tolerance = 10^(-10);
        if (abs(total-sum(ps.shunt.P)) <= tolerance)
            ResilienceTri[iterat] = 0;
        else
        ## run step 3
        Restore = RLSOPF1!(total,ps,failures,recovery_times,Pd_max);#,load_cost) # data frame [times, load shed in cost per hour]
        ## run step 4
        ResilienceTri[iterat] = crisp_res(Restore);
        end
    end
    case_res = DataFrame(resilience = ResilienceTri[:,1]);
    ## save data
    CSV.write("results\\$out_folder", case_res);

end
