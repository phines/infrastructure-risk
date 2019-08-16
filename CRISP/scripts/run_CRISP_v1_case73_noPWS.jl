#include("..\\src\\CRISP_Rdist.jl")
include("..\\src\\CRISP_Rdist.jl")
## folder of case data
case1 = "data\\saved_ps\\case73_noPWS_n-1\\"
case2 = "data\\saved_ps\\case73_noPWS\\"
out = "\\100"
if isdir("results\\"*out)
else
    mkdir("results\\"*out)
end
# number of events
N = 100000;
#set randomized seed
rng = MersenneTwister(1000);
#save restoration data to folder within results folder:
filename1 = "res_out_case39_n-1_p1"
out_folder1 = out*"\\$filename1.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist(N,case1,out_folder1)
#set randomized seed
rng = MersenneTwister(1000);
#save restoration data to folder within results folder:
filename2 = "res_out_case39_n-1_05PV_p1";
out_folder2 = out*"\\$filename2.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist(N,case2,out_folder2)
