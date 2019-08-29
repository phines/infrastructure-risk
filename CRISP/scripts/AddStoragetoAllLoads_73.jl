using CSV
using Random
include("..\\src\\CRISP_network_gen.jl")
include("..\\src\\CRISP_LSOPF_gen.jl")
include("..\\src\\CRISP_RLSOPF_gen.jl")
## name base for new case
new_case = "data\\saved_ps\\case73_noPWS_lx2_n-1+S";
## load the case data
old_case = "data\\saved_ps\\case73_noPWS_lx2_n-1\\";
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\$old_case")
ps1 = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\saved_ps\\case73_n-1")
crisp_dcpf_g_s!(ps)
total = sum(ps.shunt[:P]);
storage = deepcopy(ps.storage);
nd = length(ps.shunt[:P]);
percent_batP = 0.05;
#LoadCancList = zeros(nd);
refbus = ps.bus.id[ps.bus.bus_type.==3];
#add generators
for s = 1:nd #N_panels
    #Params = randperm(length(gen.Pg))[1];
    Bat = deepcopy(ps1.storage);
    Bat.bus .= Int64(ps.shunt.bus[s]);
    Bat.Psmax .= percent_batP.*ps.shunt.P[s]; #could also just reduce the load...
    Bat.Psmin .= -percent_batP.*ps.shunt.P[s]; #charging
    Bat.Emax .= 3 .*percent_batP.*ps.shunt.P[s]; #ratio to Power rating is based on tesla power wall
    Bat.E .= Bat.Emax;
    Bat.Emin .= 0;
    #LoadCancList[s] = percent_solar.*ps.shunt.P[s];
    Bat.Ps .= 0.0; #starts out not providing power, so model can solve the DC power flow
    println(deepcopy(ps1.storage))
    append!(storage,Bat)
end
ps.storage = storage;

PB = Int64.(100*percent_batP)
# save ps structure
if isdir(new_case*"$PB")
else
    mkdir(new_case*"$PB")
end
export_ps(ps,new_case*"$PB")
#save case info
line1 = "Original_case = CRISP\\$old_case ."
line2 = "Added battery to cover $PB fraction of demand when discharging.  "
line3 = "To each PQ buses.  "

write(new_case*"$PB\\case_info.txt",line1*line2*line3)

#include("..\\src\\CRISP_Rdist.jl")
#res = Res_dist(1000,"data\\saved_ps\\case39_05PV\\","results\\case39_05PV\\resilience_costs.csv")
#save restoration picture to:
#filename = "case39_05_gen_on_some_loads";
# save a png
#png(P,"results\\$filename")
