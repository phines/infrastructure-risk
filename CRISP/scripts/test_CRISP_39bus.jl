using CSV
include("..\\src\\CRISP_initiate.jl")
include("..\\src\\CRISP_LSOPF.jl")
include("..\\src\\CRISP_RLSOPF.jl")
include("..\\src\\CRISP_RT.jl")
include("..\\src\\CRISP_network.jl")

## number of failure scenarios to run through
Num = 4;

## load the case data
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case39\\") #case39\\")
#ps = import_ps("../data/case6ww/")
crisp_dcpf!(ps)
total = sum(ps.shunt[:P]);
Pd_max = deepcopy(ps.shunt[:P]);

# parameters of distributions for line outages and recovery times
#lines_dist = CSV.read("line-distribution-parameters.csv");
s_line = 2.56;#lines_dist[1];
maxLinesOut = length(ps.branch.f); # => k in zipf distribution
mu_line = 3.66;#lines_dist[2];
sigma_line = 2.43;#lines_dist[3];

for iterat in 1:Num
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
        crisp_dcpf!(psi);
        # run lsopf
        crisp_lsopf!(psi);
        add_changes!(ps,psi,ps_islands[i]);
        crisp_dcpf!(psi);
    end
    ## save initial outage ps file
    # make export routine, Pkg.JLD
    ## run step 3
    Restore = RLSOPF!(total,ps,failures,recovery_times,Pd_max);#,load_cost) # data frame [times, load shed in cost per hour]

    ## run step 4
    ResilienceTri = crisp_res(Restore);

    ## save data
    using CSV
    CSV.write("results\\case39\\test_initial_outage_case39_2$iterat.csv", Lines_Init_State);
    CSV.write("results\\case39\\test_restoration_case39_2$iterat.csv", Restore);
    ## make figure
    using Plots; using StatsPlots
    plot1 = @df Restore plot(:time, :perc_load_served,
            title = "Resilience Triangle",
            xlabel = "time", ylabel = "load served (%)")
    plot2 = @df Restore plot(:time, :num_lines_out,
            title = "Line Restoration",
            xlabel = "time", ylabel = "number of lines out")
    #putting 2 plots together
    P = plot(plot1,plot2,layout = (2,1),legend=false,grid=false)
    # save a png
    png(P,"results\\case39\\ResTri_case39_2$iterat")
end
