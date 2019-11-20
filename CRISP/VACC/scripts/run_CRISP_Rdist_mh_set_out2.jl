include("../src/CRISP_Rdist_vacc.jl")
#include("../src/CRISP_Rdist_PSCC_casc.jl")
## folder of case data
events = "data/outage_data/casc2/out_case73_noPWS_lx2_n-1"
case = "data/saved_ps/case73_noPWS_lx2_n-1+S5"
out = "/experiments/mh/casc2/case73_noPWS_lx2_n-1+S5"
if isdir("results"*out)
else
    mkdir("results"*out)
end
# number of events
N = 10;
# input rng iterate
i = parse(Int,ARGS[1]);
y = parse(Int,ARGS[2]);
iter = (N*100*y)+(N*i)
#time steps
dt = 60 #minutes
for j in 1:N
	#save restoration data to folder within results folder:
	filename = "res_out_$(case[20:end])";
	out_folder = out*"/$filename-$(iter+j).csv"
	# run to save csv of resilience cost distribution to the specified out_folder
	res = Resilience((iter+j),case,out_folder,events,dt)
end
