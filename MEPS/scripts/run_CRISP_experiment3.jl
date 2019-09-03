#include("..\\src\\CRISP_Rdist.jl")
include("..\\src\\CRISP_Rdist.jl")
include("..\\src\\CRISP_Rdist_no_0_OE.jl")
## folder of case data
case1 = "data\\case39\\"
# number of events
N = 1000;
#set randomized seed
rng = MersenneTwister(1000);
#save restoration data to folder within results folder:
filename1 = "res_out_case39_0_out";
out_folder1 = "\\experiments\\2\\$filename1.csv"
# run to save csv of resilience cost distribution to the specified out_folder
Res_dist(N,case1,out_folder1)
#set randomized seed
rng = MersenneTwister(0);
#save restoration data to folder within results folder:
filename2 = "res_out_case39_fair_sample_no_0_out";
out_folder2 = "\\experiments\\2\\$filename3.csv"
# run to save csv of resilience cost distribution to the specified out_folder
Res_dist_test(N,case1,out_folder2)
