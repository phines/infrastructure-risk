using CSV
using Random
include("..\\src\\CRISP_network.jl")
include("..\\src\\CRISP_LSOPF.jl")
## load the case data
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case39\\")
crisp_dcpf!(ps)
total = sum(ps.shunt[:P]);
Pd = deepcopy(ps.shunt[:P]);
gen = deepcopy(ps.gen);

#save restoration picture to:
filename = "case39_30_gen_on_some_loads";

#changing case by adding a small amount of generation at every bus with load (to be somewhat similar to roof top solar)
percent_solar = .30; #percent of demand that is served by solar for this spread of demand
perc_load_w_s = .60; #percent of load nodes that have solar on them

#pick buses with new solar
N_panels = Int64(round(perc_load_w_s*length(Pd)));
D_bus = randperm(length(ps.shunt.bus))[1:N_panels];

#pick amounts of solar
SolarP = percent_solar*sum(Pd)*rand(N_panels);

#add generators
for s = 1:N_panels
    Params = randperm(length(gen.Pg))[1];
    PV = deepcopy(gen[Params,:]);
    PV.bus = D_bus[s];
    PV.Pmax = SolarP[s];
    append!(gen,PV)
end

ps.gen = gen;

# save ps structure
export_ps(ps,"data\\saved_ps\\case39_30PV")

include("..\\src\\CRISP_Rdist.jl")
res = Res_dist(1000,"data\\saved_ps\\case39_30PV\\","results\\case39_30PV\\")
# save a png
#png(P,"results\\$filename")
