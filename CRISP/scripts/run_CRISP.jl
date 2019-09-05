include("..\\src\\CRISP_Rdist.jl")
## folder of case data
case1 = "data\\case39\\"
case2 = "data\\saved_ps\\case39+PV5\\"
case3 = "data\\saved_ps\\case39+PV20\\"
case4 = "data\\saved_ps\\case39_n-1\\"
case5 = "data\\saved_ps\\case39_n-1+PV5\\"
case6 = "data\\saved_ps\\case39_n-1+PV20\\"

out = ""
if isdir("results\\"*out)
else
    mkdir("results\\"*out)
end
# number of events
N = 10;
#set randomized seed
rng = MersenneTwister(1000);
#save restoration data to folder within results folder:
filename1 = "res_out_case39_n-1"
out_folder1 = out*"\\$filename1.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist(N,case1,out_folder1)
#set randomized seed
rng = MersenneTwister(1000);
#save restoration data to folder within results folder:
filename2 = "res_out_case39_n-1_05PV";
out_folder2 = out*"\\$filename2.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist(N,case2,out_folder2)
#set randomized seed
rng = MersenneTwister(1000);
#save restoration data to folder within results folder:
filename3 = "res_out_case39_n-1_20PV";
out_folder3 = out*"\\$filename3.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist(N,case3,out_folder3)
#set randomized seed
rng = MersenneTwister(1000);
#save restoration data to folder within results folder:
filename4 = "res_out_case39_n-1"
out_folder4 = out*"\\$filename4.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist(N,case4,out_folder4)
#set randomized seed
rng = MersenneTwister(1000);
#save restoration data to folder within results folder:
filename5 = "res_out_case39_n-1_05PV";
out_folder5 = out*"\\$filename5.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist(N,case5,out_folder5)
#set randomized seed
rng = MersenneTwister(1000);
#save restoration data to folder within results folder:
filename6 = "res_out_case39_n-1_20PV";
out_folder6 = out*"\\$filename6.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist(N,case6,out_folder6)
