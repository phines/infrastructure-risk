using CSV
using Glob
using DataFrames
#=
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
for n in 1:N
   if sum(occursin.("-$n.csv",data)) == 0
	missed[n] = true
   end
end
Index = 1:N
Index = Index[missed]
Missed = DataFrame(Missed = Index);
CSV.write("results/experiments/mh/casc2/MissedEvents_"*case*".csv",Missed)
end
=#