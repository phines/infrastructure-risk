include("..\\src\\CRISP_Rdist.jl")
## folder of case data
case = "data\\case39\\"
# number of events
N = 1000;
#set randomized seed
rng = MersenneTwister(0);
#save restoration data to folder within results folder:
filename = "resilience_out_test1";
out_folder = "\\case39\\$filename.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist(N,case,out_folder)
#set randomized seed
rng = MersenneTwister(0);
#save restoration data to folder within results folder:
filename2 = "resilience_out_test2";
out_folder2 = "\\case39\\$filename2.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist(N,case,out_folder2)
