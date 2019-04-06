#include("..\\src\\CRISP_Rdist.jl")
include("..\\src\\CRISP_Rdist_test.jl")
## folder of case data
case1 = "data\\case39\\"
case2 = "data\\saved_ps\\case39_05PV\\"
case3 = "data\\saved_ps\\case39_30PV\\"
# number of events
N = 1000;
#set randomized seed
rng = MersenneTwister(0);
#save restoration data to folder within results folder:
filename1 = "res_out_case39";
out_folder1 = "\\experiments\\1\\$filename1.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist_test(N,case1,out_folder1)
#set randomized seed
rng = MersenneTwister(0);
#save restoration data to folder within results folder:
filename2 = "res_out_case39_05PV";
out_folder2 = "\\experiments\\1\\$filename2.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist_test(N,case2,out_folder2)
#set randomized seed
rng = MersenneTwister(0);
#save restoration data to folder within results folder:
filename3 = "res_out_case39_30PV";
out_folder3 = "\\experiments\\1\\$filename3.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist_test(N,case3,out_folder3)

using CSV
using DataFrames

out1 = CSV.read("results\\"*out_folder1)
out2 = CSV.read("results\\"*out_folder2)
tol = 0.01;
if sum(abs(out1[:,1]-out2[:,1]))<=tol
println("disagree at cost on lines")
disagree = collect(1:N)[out1[:,1].!=out2[:,1]]
println(disagree[1]);
else
    println("results agree")
end
