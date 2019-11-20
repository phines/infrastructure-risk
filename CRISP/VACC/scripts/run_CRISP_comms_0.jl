include("../src/CRISP_Rdist_vacc.jl")
#include("../src/CRISP_Rdist_PSCC_comms.jl")
## folder of case data
case = "data/saved_ps/case73_noPWS_lx2_n-1+S50"
events = "data/outage_data/communication_factor/out_case73_noPWS_lx2_n-1"
out = "/experiments/mh/communications_interactions/"*case[15:end]

if isdir("results"*out)
else
    mkdir("results"*out)
end
# number of events
N = 9;
# input rng iterate
i = parse(Int,ARGS[1]);
y = parse(Int,ARGS[2]);
iter = ((N+1)*100*y)+((N+1)*i)
#time steps
dt = 60 #minutes
for i in 0:N
	#save restoration data to folder within results folder:
	filename = "res_out_$(case[20:end])";
	out_folder = out*"/$filename-$(iter+i).csv"
	# run to save csv of resilience cost distribution to the specified out_folder
	res = Resilience((iter+i),case,out_folder,events,dt)
end
