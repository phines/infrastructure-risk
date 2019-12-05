using Glob
using CSV
using DataFrames
fold = "casc2/"
fold2 = "casc\\"
l = length(fold);
Files = glob("../VACC/results/experiments/mh/"*fold*"res_*")
#include("src/CRISP_Rdist_PSCC_comms.jl")
#include("src/CRISP_Rdist_PSCC_casc.jl")
include("../src/CRISP_Rdist.jl")
#events = "data/outage_data/out_case73_noPWS_lx2_n-1"
#events = "data/outage_data/communication_factor/out_case73_noPWS_lx2_n-1"
events = "data/casc2/out_case73_noPWS_lx2_n-1"
for path in Files
	## folder of case data
	case = "data\\saved_ps\\"*path[(33+l):end-4]
	out = "\\experiments\\"*fold*path[(33+l):end-4]
	#time steps
	dt = 60 #minutes
	for m in 872
		#save restoration data to folder within results folder:
		filename = "res_out_"*path[(34+l):end-4];
		out_folder = out*"/$filename-$m.csv"
		# run to save csv of resilience cost distribution to the specified out_folder
		res = Resilience(m,case,out_folder,events,dt)
	end
end
