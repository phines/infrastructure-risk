#include("..\\src\\CRISP_Rdist.jl")
include("..\\src\\CRISP_Rdist_gen.jl")
## folder of case data
case1 = "data\\case39\\"
case2 = "data\\saved_ps\\case39+PV5\\"
case3 = "data\\saved_ps\\case39+PV20\\"
case4 = "data\\saved_ps\\case39+PV100\\"

out = "\\experiments_gen\\1"
if isdir("results\\"*out)
else
    mkdir("results\\"*out)
end

# number of events
N = 10;
#set randomized seed
rng = MersenneTwister(100);
#save restoration data to folder within results folder:
filename1 = "res_out_case39_A0O_3";
out_folder1 = out*"\\$filename1.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist_gen(N,case1,out_folder1)
#set randomized seed
rng = MersenneTwister(100);
#save restoration data to folder within results folder:
filename2 = "res_out_case39_05PV_A0O_3";
out_folder2 = out*"\\$filename2.csv"
# run to save csv of resilience cost distribution to the specified out_folder
#res = Res_dist(N,case2,out_folder2)
#set randomized seed
rng = MersenneTwister(100);
#save restoration data to folder within results folder:
filename3 = "res_out_case39_20PV_A0O_3";
out_folder3 = out*"\\$filename3.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist_gen(N,case3,out_folder3)
#set randomized seed
rng = MersenneTwister(100);
#save restoration data to folder within results folder:
filename4 = "res_out_case39_100PV_A0O_3";
out_folder4 = out*"\\$filename4.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist_gen(N,case4,out_folder4)
#using CSV
#using DataFrames

#=
out1 = CSV.read("results\\"*out_folder1, allowmissing=:none)
out2 = CSV.read("results\\"*out*"\\res_out_case39_A0O_2.csv", allowmissing=:none)
tol = 0.0001;
if sum(abs.(out1[:,1]-out2[:,1]))>tol
    println("disagree at cost on lines")
   disagree = collect(1:N)[out1[:,1].!=out2[:,1]]
   println(disagree[1]);
else
    println("results agree")
end
=#
