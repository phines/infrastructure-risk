using CSV
using Glob
using DataFrames

data = glob("results/experiments/mh/case73_n-1_l1.5/1/*0.csv");
data = [data; glob("results/experiments/mh/case73_n-1_l1.5/1/*1.csv")];
data = [data; glob("CRISP/results/experiments/mh/case73_n-1_l1.5/1/*2.csv")];
data = [data; glob("CRISP/results/experiments/mh/case73_n-1_l1.5/1/*3.csv")];
data = [data; glob("CRISP/results/experiments/mh/case73_n-1_l1.5/1/*4.csv")];
data = [data; glob("CRISP/results/experiments/mh/case73_n-1_l1.5/1/*5.csv")];
data = [data; glob("CRISP/results/experiments/mh/case73_n-1_l1.5/1/*6.csv")];
data = [data; glob("CRISP/results/experiments/mh/case73_n-1_l1.5/1/*7.csv")];
data = [data; glob("CRISP/results/experiments/mh/case73_n-1_l1.5/1/*8.csv")];
data = [data; glob("CRISP/results/experiments/mh/case73_n-1_l1.5/1/*9.csv")];

N = length(data)
r = zeros(N); 
for n in 1:N
     df = data[n]
     date = CSV.File(df) |> DataFrame
     r[n] = date.resilience[1];
end
Res = DataFrame(resilience = r);
CSV.write("res_case73_n-1_l15_1.csv",Res)