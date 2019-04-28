using CSV
include("..\\src\\CRISP_initiate.jl")
include("..\\src\\CRISP_LSOPF.jl")
include("..\\src\\CRISP_RLSOPF.jl")
include("..\\src\\CRISP_RT.jl")
include("..\\src\\CRISP_network.jl")
#set randomized seed
rng = MersenneTwister(0);

## number of failure scenarios to run through
Num = 1000;
# initialize vector of costs from events
ResilienceTri = Array{Float64}(undef,Num,1);
## load the case data
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case39\\") #case39\\")
#ps = import_ps("../data/case6ww/")
crisp_dcpf!(ps)
total = sum(ps.shunt[:P]);
Pd_max = deepcopy(ps.shunt[:P]);
ps0 = deepcopy(ps);

# parameters of distributions for line outages and recovery times
#lines_dist = CSV.read("line-distribution-parameters.csv");
s_line = 2.56;#lines_dist[1];
maxLinesOut = length(ps.branch.f); # => k in zipf distribution
mu_line = 3.66;#lines_dist[2];
sigma_line = 2.43;#lines_dist[3];

for iterat in 1:Num
    ps = deepcopy(ps0); # reset network to original state
    # step 1
    Lines_Init_State = line_state2!(ps,s_line,maxLinesOut,mu_line,sigma_line)
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
        crisp_dcpf!(psi);
        # run lsopf
        crisp_lsopf!(psi);
        ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
        ps.shunt[ps_islands[i].shunt,:P] = psi.shunt.P
        crisp_dcpf!(psi);
    end
    ## run step 3
    Restore = RLSOPF!(total,ps,failures,recovery_times,Pd_max);#,load_cost) # data frame [times, load shed in cost per hour]

    ## run step 4
    ResilienceTri[iterat] = crisp_res(Restore);
end
    case39_res = DataFrame(resilience = ResilienceTri[:,1]);
    ## save data
    CSV.write("results\\case39\\resilience_case39_5.csv", case39_res);
