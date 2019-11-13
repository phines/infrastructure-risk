using Glob
using CSV
using DataFrames
fold = "casc2/"
l = length(fold);
Files = glob("results/experiments/mh/"*fold*"MissedEvents_*")
#include("src/CRISP_Rdist_PSCC_comms.jl")
include("src/CRISP_Rdist_PSCC_casc.jl")
events = "data/outage_data/communication_factor/out_case73_noPWS_lx2_n-1"
for path in Files
	missed = CSV.File(path) |> DataFrame
	## folder of case data
	case = "data/saved_ps/"*path[(37+l):end-4]
	out = "/experiments/mh/"*fold*path[(37+l):end-4]
	#time steps
	dt = 60 #minutes
	for m in missed.Missed
		#save restoration data to folder within results folder:
		filename = "res_out_"*path[(38+l):end-4];
		out_folder = out*"/$filename-$m.csv"
		# run to save csv of resilience cost distribution to the specified out_folder
		#res = Resilience(m,case,out_folder,events,dt)
		res = Resilience(m,case,out_folder,dt)
	end
end