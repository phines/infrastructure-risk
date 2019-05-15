using CSV
using Random
include("..\\src\\CRISP_network.jl")
include("..\\src\\CRISP_LSOPF_tests.jl")
include("..\\src\\CRISP_RLSOPF_test.jl")

## load the case data
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case6ww\\")
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
    Params = randperm(length(gen.Pg))[1];
    PV = deepcopy(gen[Params,:]);
    PV.bus = D_bus[s];
    PV.Pmax = SolarP[s]; #could also just reduce the load...
    PV.Pg = 0; #starts out not providing power, so model can solve the DC power flow
    append!(gen,PV)
end

ps.gen = gen;
ps.branch.status = zeros(length(ps.branch.f));

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

ps.branch.status = 1.0.*ones(length(ps.branch.f));
subgraph = find_subgraphs(ps);# add Int64 here hide info here
M = Int64(findmax(subgraph)[1]);
ps_islands = build_islands(subgraph,ps);# at some point check for changes in islands and don't run power flows if no change
## for every island that changed (eventually)
for j in 1:M
    psi = ps_subset(ps,ps_islands[j]);
        crisp_dcpf1!(psi);
        # run lsopf
        (dPd, dPg) = crisp_rlopf1(psi,Pd[ps_islands[j].shunt]);
        # apply the results
        ps.gen[ps_islands[j].gen,:Pg]  += dPg;
        ps.shunt[ps_islands[j].shunt,:P] += dPd;
        crisp_dcpf1!(psi);
#            end
end

# save ps structure
export_ps(ps,"data\\saved_ps\\case6ww_05PV")
#save case info
line1 = "Original_case = CRISP\\data\\case6ww\\ ."
line2 = "Added solar to cover $percent_solar percent of demand.  "
line3 = "To $perc_load_w_s fraction of demand buses.  "

write("data\\saved_ps\\case6ww_05PV\\case_info.txt",line1*line2*line3)

#include("..\\src\\CRISP_Rdist.jl")
#res = Res_dist(1000,"data\\saved_ps\\case39_05PV\\","results\\case39_05PV\\resilience_costs.csv")
#save restoration picture to:
#filename = "case39_05_gen_on_some_loads";
# save a png
#png(P,"results\\$filename")
