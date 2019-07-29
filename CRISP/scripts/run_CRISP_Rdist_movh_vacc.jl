include("..\\src\\CRISP_network_gen.jl")
include("..\\src\\CRISP_Rdist_gen_movh_vacc.jl")
## folder of case data
case = "data\\saved_ps\\case39_n-1_gen\\"
case1 = "data\\saved_ps\\case39_n-1_gen+S5\\"
case2 = "data\\saved_ps\\case39_n-1_gen+S20\\"
case3 = "data\\saved_ps\\case39_n-1_gen+S50\\"

out = "\\experiments_gen_stor\\case39\\1"
if isdir("results\\"*out)
else
    mkdir("results\\"*out)
end
# iteration number
i = ARGS[1];
#time steps
dt = 10 #minutes
#set randomized seed
rng = MersenneTwister(100+$i);
#save restoration data to folder within results folder:
filename = "res_out_$(case[15:end-1])";
out_folder = out*"\\$filename$i.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist(case,out_folder,dt)


#=
#set randomized seed
rng = MersenneTwister(100);
#save restoration data to folder within results folder:
filename4 = "res_out_case73_gendata";
out_folder4 = out*"\\$filename4.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist_gen(N,case4,out_folder4)
#using CSV
#using DataFrames
=#
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
