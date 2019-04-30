#include("..\\src\\CRISP_Rdist.jl")
include("..\\src\\CRISP_Rdist.jl")
include("..\\src\\CRISP_Rdist_no_0_OE.jl")
## folder of case data
case1 = "data\\case6ww\\"
# number of events
N = 10000;
#set randomized seed
rng = MersenneTwister(1000);
#save restoration data to folder within results folder:
filename0 = "res_out_case6ww_0_out_p1";
out_folder0 = "\\experiments\\0\\$filename0.csv"
# run to save csv of resilience cost distribution to the specified out_folder
Res_dist_test2(N,case1,out_folder0)
#set randomized seed
rng = MersenneTwister(1000);
#save restoration data to folder within results folder:
filename1 = "res_out_case6ww_0_out2";
out_folder1 = "\\experiments\\0\\$filename1.csv"
# run to save csv of resilience cost distribution to the specified out_folder
Res_dist(N,case1,out_folder1)
#set randomized seed
rng = MersenneTwister(1000);
#save restoration data to folder within results folder:
filename2 = "res_out_case6ww_fair_sample_no_0_out2";
out_folder2 = "\\experiments\\0\\$filename2.csv"
# run to save csv of resilience cost distribution to the specified out_folder
Res_dist_test(N,case1,out_folder2)
