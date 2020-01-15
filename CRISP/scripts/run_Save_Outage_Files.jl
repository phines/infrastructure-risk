using Random
include("..\\src\\CRISP_Save_Outage_Files.jl")
rng = MersenneTwister(1000);
#Outages(1000,"data\\saved_ps\\case73_noPWS_lx2_n-1")
#out = "data\\outage_data\\casc_39bus\\out_case39_n-1"
out1 = "data\\outage_data\\casc_2736sp_relaxedQ"
if isdir(out1)
else mkdir(out1) end
#Outages(2,"data\\saved_ps\\case39_n-1_gen",out)
Outages(1,"data\\saved_ps\\case2736sp_relaxedQ_ps",out1*"/out_")
