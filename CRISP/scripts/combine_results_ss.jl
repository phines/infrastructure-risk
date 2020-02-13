using CSV
using Glob
using DataFrames

#include("src/CRISP_initiate.jl")
#include("src/CRISP_network.jl")


include("../src/CRISP_initiate.jl")
include("../src/CRISP_network.jl")

case = "data/saved_ps/case73_noPWS_lx2_n-1"
ps = import_ps(case)
TotalL = length(ps.branch.f)
TotalG = length(ps.gen.bus)
events = "data/strat_sample/casc_73bus" #"../data/strat_sample/casc_73bus"
Line = CSV.File(events*"/line_bins.csv") |> DataFrame
Gen = CSV.File(events*"/gen_bins.csv") |> DataFrame
L = length(Line.line_bins)
G = length(Gen.gen_bins)
Line_prob = zeros(L-1)
Gen_prob = zeros(G-1)
for l in 1:(L-1)
    bin_l = Line.line_bins[l]
    bin_h = Line.line_bins[l+1]
    Line_prob[l] = find_bins_pb(bin_l, bin_h,TotalL,true,false)
end
for g in 1:(G-1)
    bin_l = Gen.gen_bins[g]
    bin_h = Gen.gen_bins[g+1]
    Gen_prob[g] = find_bins_pb(bin_l, bin_h,TotalG,false,true)
end
load_shed_tolerance = 10^(-4);
F = glob("results/experiments/stratified_samp/*")
for fold in F # # = "results/experiments/mh/casc2+comm+ca/"
l = length(fold);
folders = glob(fold*"/case*")
#for path in folders
SS_CResults = DataFrame
path = folders[2]
case = path[l+2:end]
println(fold)
println(case)
for l in 1:(L-1)
for j in 1:(G-1)
n = Line.line_bins[l]
g = Gen.gen_bins[j]
println("$n-$g")
data = glob(path*"/$n-$g/*-*0.csv");
data = [data; glob(path*"/$n-$g/*-*1.csv")];
data = [data; glob(path*"/$n-$g/*-*2.csv")];
data = [data; glob(path*"/$n-$g/*-*3.csv")];
data = [data; glob(path*"/$n-$g/*-*4.csv")];
data = [data; glob(path*"/$n-$g/*-*5.csv")];
data = [data; glob(path*"/$n-$g/*-*6.csv")];
data = [data; glob(path*"/$n-$g/*-*7.csv")];
data = [data; glob(path*"/$n-$g/*-*8.csv")];
data = [data; glob(path*"/$n-$g/*-*9.csv")];
N = length(data) #100 # expected length of data if no files missing
missed = falses(N);
r = zeros(N);
sizes = zeros(N);
restore_time = zeros(N);
l_bins = n .*ones(N);
g_bins = g .*ones(N);
m=0;
for df in data
   m=m+1;
   #df = data[occursin.("-$m.csv",data)]
   df2 = df[1:end-4]*"_restore.csv"
   Restore = CSV.File(df2) |> DataFrame
   sizes[m] = Restore.load_shed[2];
   # find load shed larger than load_shed_tolerance
   s = Restore.load_shed .>= load_shed_tolerance
   if sum(s) .!=0
      # latest time with load shed of more than load_shed_tolerance
      restore_time[m] = Restore.time[s][end]
   end
   T = (Restore.time[3:end]-Restore.time[2:end-1])./60;
   r[m] = sum(T.*Restore.load_shed[2:end-1]);
end
Index = 1:N
Index = Index[missed];
Results = DataFrame(Resilience = r, Size = sizes, Restoration_time = restore_time, lines_bin = n .*ones(N), gens_bin = g .*ones(N), event_num = 1:N);
sort!(Results)
P = (N:-1:1)/N;
P = (P.* Line_prob[l] .* Gen_prob[j])
Results[!, :Prob] = P
if ((l==1) & (g==1))	
	SS_CResults = Results
else	
	append!(SS_CResults,Results)	
end
CSV.write(fold*"/RiskResults_"*case*"-$n-$g.csv",Results)
end
end
sort!(SS_CResults)
CSV.write(fold*"/0_CombinedResults_"*case*".csv",SS_CResults)
end
#end
#end

