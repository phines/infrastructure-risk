#include("..\\src\\CRISP_Rdist.jl")
include("..\\src\\CRISP_Rdist.jl")
## folder of case data
case1 = "data\\saved_ps\\case73_noPWS_n-1\\"
case2 = "data\\saved_ps\\case73_noPWS\\"
out = "100\\case73_load1.5"
if isdir("results\\"*out)
else
    mkdir("results\\"*out)
end
# number of events
N = 10000;
#set randomized seed
rng = MersenneTwister(1000);
#save restoration data to folder within results folder:
filename1 = "res_out_case73_n-1_p1"
out_folder1 = out*"\\$filename1.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist(N,case1,out_folder1)
