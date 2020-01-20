using Random
include("..\\src\\CRISP_initiate.jl")
rng = MersenneTwister(100);
#=
#out1 = "data/outage_data/casc_73_gt"
out1 = "data\\outage_data\\casc_39bus_gt"
#out1 = "data\\outage_data\\casc_2736sp_relaxedQ"
if isdir(out1)
else mkdir(out1) end
#case = "data\\saved_ps\\case73_noPWS_lx2_n-1"
case = "data\\saved_ps\\case39_n-1_gen"
#out = out1*"\\out_case73_noPWS_lx2_n-1"
out = out1*"\\out_case39_n-1"
Outages(1000,case,out,true)

#Polish case
#case = "data\\saved_ps\\case2736sp_relaxedQ_ps"
#Outages(1,case,out1*"/out_",false)
#Outages(1,case,out1*"/out_",true)
=#

## Stratified sampling
out1 = "data\\strat_sample\\casc_39bus_gt"
#out1 = "data\\outage_data\\casc_2736sp_relaxedQ"
if isdir(out1)
else mkdir(out1) end
#case = "data\\saved_ps\\case73_noPWS_lx2_n-1"
case = "data\\saved_ps\\case39_n-1_gen"
#out = out1*"\\out_case73_noPWS_lx2_n-1"
outfile = "/out_case39_n-1"
Outages_ss(100,case,out1,outfile,false)
