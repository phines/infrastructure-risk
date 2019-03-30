#set up packages
using CSV; using DataFrames; using SpecialFunctions;using SparseArrays
rng = MersenneTwister(0);
# uses the distributions of the number of lines outages and the distribution of recovery times
# lines which are fits to BPA data to realize a specifc number of line failures and the recovery
# time of each failure
# it also uses an exponential function with paamter lambda=1 to determine the number of generators
# that are outaged
include("..\\src\\CRISP_network.jl")
include("..\\src\\CRISP_initiate_v2.jl")
# network case to model
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case39\\")
# parameters of distributions for line outages and recovery times
#lines_dist = CSV.read("line-distribution-parameters.csv");
s_line = 2.56;#lines_dist[1];
mu_line = 3.66;#lines_dist[2];
sigma_line = 2.43;#lines_dist[3];
# number of events to realize the outages and restoration times for
Nevents = 1000;
# number of lines and generators in network case
TotalLines = length(ps.branch[1]);
TotalGens =length(ps.gen[1]);

maxLinesOut = TotalLines;
Nlines = zeros(Nevents);
Times = [];
for n = 1:Nevents
    Nlines[n] = init_out_zipf(s_line,maxLinesOut,TotalLines);
end
for n = 1:Nevents
    push!(Times,RecoveryTimes(mu_line,sigma_line,1)[1,1]);
end
Nline = sort(Nlines);
ccdf_nlines_empir = collect(length(Nline):-1:1)./length(Nline);
Times = Times[:];
RecT = sort(Times);
ccdf_rectime_empir = collect(length(RecT):-1:1)./length(RecT);

num = DataFrame(Realizations_of_NumLines = Nlines, Number_Lines_Out = Nline, CCDF = ccdf_nlines_empir);
CSV.write("results\\checks\\LineOutageRealiz_v2.csv",num);
rec = DataFrame(Realizations_of_RecovTime = Times[:], Recovery_Time = RecT, CCDF = ccdf_rectime_empir);
CSV.write("results\\checks\\RecTimeRealiz_v2.csv",rec);
k = 1:100;
pdf_nlines = (Nevents/zeta(s_line)).*k.^(-s_line);
H_k_s = zeros(k[end]);
for i = 1:k[end]
    for j = 1:i
    H_k_s[i] = H_k_s[i] + 1/(j^s_line);
    end
end
ccdf_lines = (1 .- H_k_s./zeta(s_line));
NumLinesAnalytic = DataFrame(k=k, pdf_Num_Lines_Out=pdf_nlines, ccdf_Num_Lines_Out=ccdf_lines);
CSV.write("results\\checks\\LinesDist_v2.csv",NumLinesAnalytic);
T = 0:0.1:10000;
ccdf_rectime = zeros(length(T));
for i = 1:length(T)
    t = T[i];
    ccdf_rectime[i] = 1-(0.5(1+erf((log.(t)-mu_line)./(sqrt(2)*sigma_line))));
end
pdf_rectime = ((length(Times) .- 1).*(1 ./ (sigma_line.*sqrt.(2 .*pi)))).*T.^(-1).*exp.(-((log.(T).- mu_line).^2).*(1 ./ (2 .*sigma_line.^2)));
RecovTimeAnalytic = DataFrame(T = T, pdf_RecoveryTime=pdf_rectime, ccdf_RecoveryTime=ccdf_rectime)
CSV.write("results\\checks\\RecTimeDist_v2.csv",RecovTimeAnalytic);
#using Gadfly
#NumLines = plot(layer(NumLinesAnalytic, x="k", y="ccdf_Num_Lines_Out", Theme(default_color=colorant"red"), Geom.line), layer(num, x = "Number_Lines_Out", y ="CCDF", Geom.bar), Scale.x_log10);
#RestoreT = plot(layer(RecovTimeAnalytic, x="T", y="ccdf_RecoveryTime", Theme(default_color=colorant"red"), Geom.line), layer(rec, x = "Recovery_Time", y = "CCDF", Geom.bar), Scale.x_log10, Scale.y_log10);
#draw(PDF("validate-BPA-line-outage.pdf", 10inch, 6inch), NumLines);
#draw(PDF("validate-BPA-recovery-time.pdf", 10inch, 6inch), RestoreT);
