using Random
include("..\\src\\CRISP_initiate.jl")
rng = MersenneTwister(100);

out1 = "data/outage_data/casc_73"
#out1 = "data\\outage_data\\casc_39bus_gt"
#out1 = "data\\outage_data\\casc_2736sp_relaxedQ"
if isdir(out1)
else mkdir(out1) end
case = "data\\saved_ps\\case73_noPWS_lx2_n-1"
#case = "data\\saved_ps\\case39_n-1_gen"
out = out1*"\\out_case73_noPWS_lx2_n-1"
#out = out1*"\\out_case39_n-1"
Outages(10000,case,out)
#=
#Polish case
#case = "data\\saved_ps\\case2736sp_relaxedQ_ps"
#Outages(1,case,out1*"/out_",false)
#Outages(1,case,out1*"/out_",true)
=#

## Stratified sampling
#=
## 73 bus case
out1 = "data\\strat_sample_bins\\casc_73bus"
if isdir(out1)
else mkdir(out1) end
case = "data\\saved_ps\\case73_noPWS_lx2_n-1"
outfile = "/out_case73_noPWS_lx2_n-1"
n = 100
l_bins = 20
g_bins = 5
Outages_ss(n,case,out1,outfile,l_bins,g_bins,false)
=#
#=
##39 bus case
out1 = "data\\strat_sample_bins\\casc_39bus"
if isdir(out1)
else mkdir(out1) end
case = "data\\saved_ps\\case39_n-1_gen"
outfile = "/out_case39_n-1"
n = 100
l_bins = 20
g_bins = 5
Outages_ss(n,case,out1,outfile,l_bins,g_bins,false)
=#
#=
## Polish case
#out1 = "data\\strat_sample_bins\\casc_39bus_gt"
out1 = "data\\strat_sample_bins\\casc_73bus"
#out1 = "data\\strat_sample\\casc_2736sp_relaxedQ"
if isdir(out1)
else mkdir(out1) end
case = "data\\saved_ps\\case73_noPWS_lx2_n-1"
#case = "data\\saved_ps\\case39_n-1_gen"
#case = out_case73_noPWS_lx2_n-1"
#outfile = "/out_case39_n-1"
outfile = "/out_case73_noPWS_lx2_n-1"
n = 20
l_bins = 20
g_bins = 5
Outages_ss(n,case,out1,outfile,l_bins,g_bins,false)
=#
