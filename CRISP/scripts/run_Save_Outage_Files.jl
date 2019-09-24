using Random
include("..\\src\\CRISP_Save_Outage_Files.jl")
rng = MersenneTwister(100);
Outages(1000,"data\\saved_ps\\case73_noPWS_lx2_n-1")
