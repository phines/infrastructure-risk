include("..\\src\\CRISP_Rdist_test1.jl")
include("..\\src\\CRISP_Rdist_test.jl")
## folder of case data
case = "data\\case39\\"
# number of events
N = 1000;
#set randomized seed
rng = MersenneTwister(0);
#save restoration data to folder within results folder:
filename1 = "res_out_test_008";
out_folder1 = "\\case39\\test\\$filename1.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist_test1(N,case,out_folder1)
#set randomized seed
rng = MersenneTwister(0);
#save restoration data to folder within results folder:
filename2 = "res_out_test_009";
out_folder2 = "\\case39\\test\\$filename2.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist_test(N,case,out_folder2)

using CSV
using DataFrames

out1 = CSV.read("results\\"*out_folder1)
out2 = CSV.read("results\\"*out_folder2)
tol = 0.0001;
if sum(abs.(out1[:,1]-out2[:,1]))>tol
    println("disagree at cost on lines")
    disagree = collect(1:N)[out1[:,1].!=out2[:,1]]
    println(disagree[1]);
else
    println("results agree")
end
