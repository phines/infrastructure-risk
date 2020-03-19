using Random
include("..\\src\\CRISP_initiate.jl")
rng = MersenneTwister(100);

## Stratified sampling

## 73 bus case
out1 = "data\\strat_sample_bins\\NE_NY"
if isdir(out1)
else mkdir(out1) end
case = "../../../NE_NY_data/ps_reduced_bs"
outfile = "/out_ne_ny"
n = 100
l_bins = 20
g_bins = 5
Outages_ss(n,case,out1,outfile,l_bins,g_bins,false;bins_lines=[],bins_gens=[0 1 3 5 10 20])
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
