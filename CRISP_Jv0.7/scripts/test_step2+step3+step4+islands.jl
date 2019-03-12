#using CRISP_LSOPF
using CSV
#include code for all necessary steps (2,3,4,and grid segmenting)
include("..\\src\\CRISP_LSOPF.jl")
include("..\\src\\CRISP_RLSOPF.jl")
include("..\\src\\CRISP_RT.jl")
include("..\\src\\CRISP_network_segments.jl")
include("..\\src\\parser2.jl")
include("..\\src\\s1-initiate2.jl")
#include("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\src\\parser.jl")
## load the case data
ps = mp2ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case6ww.m"); #case39.m")
#ps = mp2ps("../data/case6ww.m")
crisp_dcpf!(ps)
total = sum(ps.shunt[:P]);
Pd_max = deepcopy(ps.shunt[:P]);

# remove branches that result in islands
Nlines = 6;
ps.branch[2,:status]=0;
ps.branch[4,:status]=0;
ps.branch[5,:status]=0;
ps.branch[7,:status]=0;
ps.branch[8,:status]=0;
ps.branch[11,:status]=0;

failures = ps.branch.status;
mu_line = 3.66;#lines_dist[2];
sigma_line = 2.43;#lines_dist[3];
RecovTimeL = RecoveryTimes(mu_line,sigma_line,Nlines);
lines_outage_recovery = RecTime(RecovTimeL,failures);
recovery_times = lines_outage_recovery[:,2];
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
## save initial outage ps file
# make export routine, Pkg.JLD
## run step 3
Restore = RLSOPF!(total,ps,failures,recovery_times,Pd_max);#,load_cost) # data frame [times, load shed in cost per hour]

## run step 4
ResilienceTri = crisp_res(Restore);

## save data
using CSV
CSV.write("results\\test_initial_outage_case6ww.csv", Lines_Init_State);
CSV.write("results\\test_restoration_case6ww.csv", Restore);
## make figure
using StatsPlots
%#df Restore plot(:time, :load_shed,
#        title = "Resilience Triangle",
#        xlabel = "time", ylabel = "load shed")

@df Restore plot(:time, :perc_load_served,
                title = "Resilience Triangle",
                xlabel = "time", ylabel = "percent load served")

# save a png
png("results\\ResTri_perc")
