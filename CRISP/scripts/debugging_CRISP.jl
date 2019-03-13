# Debugging results/case39/*case39_9.*

using CSV

include("..\\src\\CRISP_LSOPF.jl")
include("..\\src\\CRISP_RLSOPF.jl")
include("..\\src\\CRISP_RT.jl")
include("..\\src\\CRISP_network.jl")
## load the case data
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case39\\")
crisp_dcpf!(ps)
total = sum(ps.shunt[:P]);
Pd_max = deepcopy(ps.shunt[:P]);

#save restoration picture to:
filename = "Fixed_case39_9";

#pull failures from csv file
Lines_Init_State = CSV.read("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\results\\Extreme_initial_outages_case39_9.csv");
state = Lines_Init_State[:,1];
recovery_times = Array{Float64}(Lines_Init_State[1:end,2]);
failures = Array{Int64}(state[1:end]);
ps.branch.status = failures;

#check for islands
subgraph = find_subgraphs(ps);
M = Int64(findmax(subgraph)[1]);
ps_islands = build_islands(subgraph,ps);
for i in 1:M
    psi = ps_subset(ps,ps_islands[i]);
    ## run step 2
    crisp_dcpf!(psi);
    # run lsopf
    crisp_lsopf!(psi);
    add_changes!(ps,psi,ps_islands[i]);
    crisp_dcpf!(psi);
end

# save ps structure
export_ps(ps,"data\\saved_ps\\case39_severe")

## run step 3
Restore = RLSOPF!(total,ps,failures,recovery_times,Pd_max);#,load_cost) # data frame [times, load shed in cost per hour]
#make_plots(Restore,filename)
using Plots; using StatsPlots
#theme(:sand)
plot1 = @df Restore plot(:time, :perc_load_served,
      title = "Resilience Triangle",
      xlabel = "time", ylabel = "load served (%)") #, legend = off, grid = off
plot2 = @df Restore plot(:time, :num_lines_out,
            title = "Line Restoration",
            xlabel = "time", ylabel = "number of lines out")
#putting 2 plots together
P = plot(plot1,plot2,layout = (2,1),legend=false,grid=false)
# save a png
png(P,"results\\$filename")
## run step 4
ResilienceTri = crisp_res(Restore);
