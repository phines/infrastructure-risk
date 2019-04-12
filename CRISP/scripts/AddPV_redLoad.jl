using CSV
using Random
include("..\\src\\CRISP_network.jl")
include("..\\src\\CRISP_LSOPF_tests.jl")
include("..\\src\\CRISP_RLSOPF_test.jl")

## load the case data
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case39\\")
crisp_dcpf1!(ps)
total = sum(ps.shunt[:P]);
Pd = deepcopy(ps.shunt[:P]);
gen = deepcopy(ps.gen);

#changing case by adding a small amount of generation at every bus with load (to be somewhat similar to roof top solar)
percent_solar = .05; #percent of demand that is served by solar for this spread of demand
perc_load_w_s = .4; #percent of load nodes that have solar on them

#pick buses with new solar
N_panels = Int64(round(perc_load_w_s*length(Pd)));
D_bus = ps.shunt.bus[randperm(length(ps.shunt.bus))];

#pick amounts of solar
InitP = rand(N_panels);
SolarP = percent_solar*sum(Pd)*InitP/sum(InitP);

#add generators
for s = 1:N_panels
    ps.shunt.P .+= -SolarP[s]; #could also just reduce the load...
end

# save ps structure
export_ps(ps,"data\\saved_ps\\case39_05PV_RS")
#save case info
line1 = "Original_case = CRISP\\data\\case39\\ ."
line2 = "Added solar to cover $percent_solar percent of demand.  "
line3 = "To $perc_load_w_s fraction of demand buses.  "
line4 = "By reducing random loads.  "

write("data\\saved_ps\\case6ww_05PV\\case_info.txt",line1*line2*line3*line4)

#include("..\\src\\CRISP_Rdist.jl")
#res = Res_dist(1000,"data\\saved_ps\\case39_05PV\\","results\\case39_05PV\\resilience_costs.csv")
#save restoration picture to:
#filename = "case39_05_gen_on_some_loads";
# save a png
#png(P,"results\\$filename")
