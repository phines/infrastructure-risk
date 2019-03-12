#using CRISP_LSOPF
using CSV
#include code for all necessary steps (2,3,4,and grid segmenting)
include("..\\src\\CRISP_LSOPF.jl")
include("..\\src\\CRISP_RLSOPF.jl")
include("..\\src\\CRISP_RT.jl")
include("..\\src\\CRISP_network.jl")
include("..\\src\\CRISP_initiate.jl")
#include("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\src\\parser.jl")
## load the case data
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case39\\")
crisp_dcpf!(ps)
total = sum(ps.shunt[:P]);
Pd_max = deepcopy(ps.shunt[:P]);

# parameters of distributions for line outages and recovery times
#lines_dist = CSV.read("line-distribution-parameters.csv");
s_line = 2.56;#lines_dist[1];
maxLinesOut = length(ps.branch.f); # => k in zipf distribution
mu_line = 3.66;#lines_dist[2];
sigma_line = 2.43;#lines_dist[3];
# parameters of distributions for generator outages and recovery times
#gens_dist = CSV.read("gen-distribution-parameters.csv");
lambda_gen = 1;#gens_dist[1];
mu_gen = 3.66;#gens_dist[2];
sigma_gen = 2.43;#gens_dist[3];

#import CRISP
# step 1
#Lines_Init_State = CRISP.line_state(ps,s_line,maxLinesOut,mu_line,sigma_line,orignumLines)
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
    ## run step 2
    # run the dcpf
    crisp_dcpf!(psi);
    # run lsopf
    (dPd, dPg) = crisp_lsopf(psi);
    # apply the results
    ps.gen[ps_islands[i].gen,:Pg]  += dPg;
    ps.shunt[ps_islands[i].shunt,:P] += dPd;
    crisp_dcpf!(psi);
end
## run step 3
Restore = RLSOPF!(total,ps,failures,recovery_times,Pd_max);#,load_cost) # data frame [times, load shed in cost per hour]

## run step 4
ResilienceTri = crisp_res(Restore);

## save data
using CSV
CSV.write("results\\test_initial_outage_case39.csv", Lines_Init_State);
CSV.write("results\\test_restoration_case39.csv", Restore);
## make figure
using StatsPlots
@df Restore plot(:time, :load_shed,
        title = "Resilience Trianglge",
        xlabel = "time", ylabel = "load shed")

# save a png
png("results\\Restri2")
#using Gadfly
# make plots
#A = Gadfly.plot(Restore, x="time", y="load_shed", Geom.line);
#save plots
#draw(PNG("results\\ResTri1.png", 6inch, 6inch), A);
