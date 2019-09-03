include("../src/CRISP_Rdist_gen_mh_varL+PV.jl")
## folder of case data
#case = "data/saved_ps/case39_n-1_gen"
#case1 = "data/saved_ps/case39_n-1_gen+S5"
#case2 = "data/saved_ps/case39_n-1_gen+S20"
#case3 = "data/saved_ps/case39_n-1_gen+S50"

case = "data\\saved_ps\\case73_noPWS_lx2_n-1+PV5"

out = "\\experiments\\mh\\var\\case73_lx2_n-1_+PV5"
if isdir("VACC\\results"*out)
else
    mkdir("VACC\\results"*out)
end
# number of events
N = 1000;
# input rng iterate
i = 0 #parse(Int,ARGS[1]);
y = 0 #parse(Int,ARGS[2]);
#iter = (100*y)+(N*i)
#time steps
dt = 60 #minutes
#set randomized seed
rng = MersenneTwister(100);#100+iter
#save restoration data to folder within results folder:
filename = "res_out_$(case[20:end])";
out_folder = out*"\\$filename.csv"#$filename-$iter.csv"
# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist(N,case,out_folder,dt)
