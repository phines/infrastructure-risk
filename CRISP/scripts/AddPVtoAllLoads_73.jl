using CSV
using Random
include("..\\src\\CRISP_network_gen.jl")
include("..\\src\\CRISP_LSOPF_gen.jl")
include("..\\src\\CRISP_RLSOPF_gen.jl")
## name base for new case
new_case = "data\\saved_ps\\case73_noPWS_lx2_n-1+PV";
## load the case data
old_case = "data\\saved_ps\\case73_noPWS_lx2_n-1\\";
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\$old_case")
ps1 = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\saved_ps\\case73_n-1")
#make_types_general!(ps)
#make_types_general!(ps1)
crisp_dcpf_g!(ps)
total = sum(ps.shunt.P);
Pd = deepcopy(ps.shunt.P);
gen = deepcopy(ps.gen);
nd = length(ps.shunt.P);
ng = length(ps.gen.Pg);
percent_solar = 0.5;
#LoadCancList = zeros(nd);
refbus = ps.bus.id[ps.bus.bus_type.==3];
#solar_options = ps1.gen.[97:116]
#solar_options = append!(solar_options,ps1.gen.[118:97+56])
#add generators
PV = deepcopy(ps1.gen[97:116,:]);
Numbr = Int64(round(nd./20));
if Numbr .>1
    for n in 1:Numbr
        append!(PV,PV);
    end
end
for s = 1:nd #N_panels
    if ps.shunt.bus[s]==refbus[1]
        ps.shunt.P[s] -= percent_solar.*ps.shunt.P[s];
    else
        PV.bus[s] = ps.shunt.bus[s];
        PV.Pmax[s] = percent_solar.*ps.shunt.P[s]; #could also just reduce the load...
        PV.Qg[s] = 0.0;
        PV.Pmin[s] = 0.0;
        #LoadCancList[s] = percent_solar.*ps.shunt.P[s];
        PV.Pg[s] = percent_solar.*ps.shunt.P[s]; #starts out not providing power, so model can solve the DC power flow
    end
    append!(gen,PV[1:nd,:]);
end
ps.gen = gen;
println(ps.gen.Pg)
crisp_dcpf_g!(ps)
println(ps.gen.Pg)
crisp_lsopf_g!(ps)
println(ps.gen.Pg)
#check
crisp_dcpf_g!(ps)
println(ps.gen.Pg)

PS = Int64.(100*percent_solar)
# save ps structure
if isdir(new_case*"$PS")
else
    mkdir(new_case*"$PS")
end
export_ps(ps,new_case*"$PS")
#save case info
line1 = "Original_case = CRISP\\$old_case ."
line2 = "Added solar to cover $percent_solar fraction of demand.  "
line3 = "To each PQ buses.  "

write(new_case*"$PS\\case_info.txt",line1*line2*line3)

#include("..\\src\\CRISP_Rdist.jl")
#res = Res_dist(1000,"data\\saved_ps\\case39_05PV\\","results\\case39_05PV\\resilience_costs.csv")
#save restoration picture to:
#filename = "case39_05_gen_on_some_loads";
# save a png
#png(P,"results\\$filename")
function make_types_general!(ps)
    ps.gen.Pmax .= ps.gen.Pmax .* 1.0
    ps.gen.Pmin .= ps.gen.Pmin .* 1.0
    ps.gen.Pg .= ps.gen.Pg .* 1.0
    ps.gen.Qg .= ps.gen.Qg .* 1.0
    return ps
end
