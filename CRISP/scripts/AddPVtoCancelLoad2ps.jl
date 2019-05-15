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
nd = length(ps.shunt[:P]);
ng = length(ps.gen.Pg);
percent_solar = 1;
#add generators
for s = 1:nd #N_panels
    Params = randperm(length(gen.Pg))[1];
    PV = deepcopy(gen[Params,:]);
    PV.bus = ps.shunt.bus[s];
    PV.Pmax = percent_solar.*ps.shunt.P[s]; #could also just reduce the load...
    PV.Pg = ps.shunt.P[s]; #starts out not providing power, so model can solve the DC power flow
    append!(gen,PV)
end
ps.gen = gen;
for g = 1:ng
    ps.gen.Pg[g] = 0.0
end
println(ps.gen.Pg)
#check
crisp_dcpf1!(ps)
println(ps.gen.Pg)
#crisp_lsopf1!(ps)
#println(ps.gen.Pg)

# save ps structure
export_ps(ps,"data\\saved_ps\\case6ww_100PV")
#save case info
line1 = "Original_case = CRISP\\data\\case6ww\\ ."
line2 = "Added solar to cover $percent_solar fraction of demand.  "
line3 = "To each Q buses.  "

write("data\\saved_ps\\case6ww_100PV\\case_info.txt",line1*line2*line3)

#include("..\\src\\CRISP_Rdist.jl")
#res = Res_dist(1000,"data\\saved_ps\\case39_05PV\\","results\\case39_05PV\\resilience_costs.csv")
#save restoration picture to:
#filename = "case39_05_gen_on_some_loads";
# save a png
#png(P,"results\\$filename")
