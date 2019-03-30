include("..\\src\\CRISP_Rdist.jl")
## folder of case data
case = "data\\case39\\"
# number of events
N = 1000;
#save restoration data to folder within results folder:
filename = "resilience_out3_copy";
out_folder = "\\case39\\$filename.csv"

# run to save csv of resilience cost distribution to the specified out_folder
res = Res_dist(N,case,out_folder)
