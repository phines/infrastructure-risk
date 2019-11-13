using CSV
using Glob
using DataFrames
case = "case73_noPWS_lx2_n-1"
data = glob("../results/experiments/mh/var/"*case*"/*0.csv");
data = [data; glob("../results/experiments/mh/var/"*case*"/*1.csv")];
data = [data; glob("../results/experiments/mh/var/"*case*"/*2.csv")];
data = [data; glob("../results/experiments/mh/var/"*case*"/*3.csv")];
data = [data; glob("../results/experiments/mh/var/"*case*"/*4.csv")];
data = [data; glob("../results/experiments/mh/var/"*case*"/*5.csv")];
data = [data; glob("../results/experiments/mh/var/"*case*"/*6.csv")];
data = [data; glob("../results/experiments/mh/var/"*case*"/*7.csv")];
data = [data; glob("../results/experiments/mh/var/"*case*"/*8.csv")];
data = [data; glob("../results/experiments/mh/var/"*case*"/*9.csv")];

N = length(data)
r = zeros(N); 
for n in 1:N
     df = data[n]
     date = CSV.File(df) |> DataFrame
     r[n] = date.resilience[1];
end
Res = DataFrame(resilience = r);
CSV.write("results/experiments/mh/var/res_"*case*".csv",Res)

#=
folders = glob("results/experiments/mh/set/case*")
for path in folders
case = path[28:end];
data = glob("results/experiments/mh/set/"*case*"/*0.csv");
data = [data; glob(path*"/*1.csv")];
data = [data; glob(path*"/*2.csv")];
data = [data; glob(path*"/*3.csv")];
data = [data; glob(path*"/*4.csv")];
data = [data; glob(path*"/*5.csv")];
data = [data; glob(path*"/*6.csv")];
data = [data; glob(path*"/*7.csv")];
data = [data; glob(path*"/*8.csv")];
data = [data; glob(path*"/*9.csv")];
N = length(data)
r = zeros(N); 
for n in 1:N
     df = data[n]
     date = CSV.File(df) |> DataFrame
     r[n] = date.resilience[1];
end
Res = DataFrame(resilience = r);
CSV.write("results/experiments/mh/set/res_"*case*".csv",Res)
end
=#
#=
folders = glob("results/experiments/mh/communications_interactions/case*")
for path in folders
case = path[56:end];
data = glob("results/experiments/mh/communications_interactions/"*case*"/*0.csv");
data = [data; glob(path*"/*1.csv")];
data = [data; glob(path*"/*2.csv")];
data = [data; glob(path*"/*3.csv")];
data = [data; glob(path*"/*4.csv")];
data = [data; glob(path*"/*5.csv")];
data = [data; glob(path*"/*6.csv")];
data = [data; glob(path*"/*7.csv")];
data = [data; glob(path*"/*8.csv")];
data = [data; glob(path*"/*9.csv")];
N = length(data)
r = zeros(N); 
for n in 1:N
     df = data[n]
     date = CSV.File(df) |> DataFrame
     r[n] = date.resilience[1];
end
Res = DataFrame(resilience = r);
CSV.write("results/experiments/mh/communications_interactions/res_"*case*".csv",Res)
end
=#
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
N = length(data)
r = zeros(N); 
for n in 1:N
     df = data[n]
     date = CSV.File(df) |> DataFrame
     r[n] = date.resilience[1];
end
Res = DataFrame(resilience = r);
CSV.write("results/experiments/mh/casc2/res_"*case*".csv",Res)
end
=#
