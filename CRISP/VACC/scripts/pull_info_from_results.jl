using CSV
using Glob
using DataFrames
load_shed_tolerance = 10^(-4);
folders = glob("results/experiments/mh/casc2/case*")
for path in folders
case = path[30:end];
data = glob("results/experiments/mh/casc2/"*case*"/*0.csv");
data = [data; glob(path*"/*1.csv")];
data = [data; glob(path*"/*2.csv")];
data = [data; glob(path*"/*3.csv")];
data = [data; glob(path*"/*4.csv")];
data = [data; glob(path*"/*5.csv")];
data = [data; glob(path*"/*6.csv")];
data = [data; glob(path*"/*7.csv")];
data = [data; glob(path*"/*8.csv")];
data = [data; glob(path*"/*9.csv")];
N = 1000 # expected length of data if no files missing
missed = falses(N);
r = zeros(N);
sizes = zeros(N);
restore_time = zeros(N);
for n in 1:N
   df = data[occursin.("-$n.csv",data)]
   df2 = df[1][1:end-4]*"_restore.csv"
   date = CSV.File(df[1]) |> DataFrame
   r[n] = date.resilience[1];
   Restore = CSV.File(df2) |> DataFrame
   sizes[n] = Restore.load_shed[2];
   # find load shed larger than load_shed_tolerance
   s = Restore.load_shed .>= load_shed_tolerance
   if sum(s) .!=0
      # latest time with load shed of more than load_shed_tolerance
      restore_time[n] = Restore.time[s][end]
   end
end
Index = 1:N
Index = Index[missed];
Results = DataFrame(Resilience = r, Size = sizes, Restoration_time = restore_time);
CSV.write("results/experiments/mh/casc2/RiskResults_"*case*".csv",Results)
end