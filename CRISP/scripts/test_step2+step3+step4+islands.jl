using CSV
include("..\\src\\CRISP_initiate.jl")
include("..\\src\\CRISP_LSOPF.jl")
include("..\\src\\CRISP_RLSOPF.jl")
include("..\\src\\CRISP_RT.jl")
include("..\\src\\CRISP_network.jl")

## load the case data
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case6ww\\") #case39\\")
#ps = import_ps("../data/case6ww/")
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
