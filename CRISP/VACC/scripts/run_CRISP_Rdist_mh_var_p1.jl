include("../src/CRISP_Rdist_vacc.jl")
#include("../src/CRISP_Rdist_gen_mh_varL+PV.jl")
## folder of case data
#case = "data/saved_ps/case39_n-1_gen"
#case1 = "data/saved_ps/case39_n-1_gen+S5"
#case2 = "data/saved_ps/case39_n-1_gen+S20"
#case3 = "data/saved_ps/case39_n-1_gen+S50"

case = "data/saved_ps/case73_noPWS_lx2_n-1+S20"

out = "/experiments/mh/var/case73_noPWS_lx2_n-1+S20"
if isdir("results"*out)
else
    mkdir("results"*out)
end
# number of events
N = 1;
# input rng iterate
i = parse(Int,ARGS[1]);
y = parse(Int,ARGS[2]);
iter = (100*y)+(N*i)
#time steps
dt = 60 #minutes
#set randomized seed
rng = MersenneTwister(100+iter);
#save restoration data to folder within results folder:
filename = "res_out_$(case[20:end])";
out_folder = out*"/$filename-$iter.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Resilience_1(N,case,out_folder,dt)
