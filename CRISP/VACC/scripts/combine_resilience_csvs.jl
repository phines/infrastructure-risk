using CSV
using Glob
using DataFrames

data = glob("CRISP/results/experiments/mh/case73_n-1/3/*-1.csv");
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*11.csv")]
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*0.csv")];
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*01.csv")]
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*2.csv")];
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*21.csv")]
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*3.csv")];
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*31.csv")]
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*4.csv")];
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*41.csv")]
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*5.csv")];
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*51.csv")]
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*6.csv")];
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*61.csv")]
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*7.csv")];
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*71.csv")]
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*8.csv")];
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*81.csv")]
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*9.csv")];
data = [data; glob("CRISP/results/experiments/mh/case73_n-1/3/*91.csv")]

N = length(data)
r = zeros(N); 
for n in 1:N
     df = data[n]
     date = CSV.File(df) |> DataFrame
     r[n] = date.resilience[1];
end
Res = DataFrame(resilience = r);
CSV.write("CRISP/results/experiments/mh/resilience_case73_n-1_3.csv",Res)