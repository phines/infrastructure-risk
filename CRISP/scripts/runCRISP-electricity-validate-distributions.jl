#set up packages
using CSV; using DataFrames; using SpecialFunctions;

# uses the distributions of the number of lines outages and the distribution of recovery times
# lines which are fits to BPA data to realize a specifc number of line failures and the recovery
# time of each failure
# it also uses an exponential function with paamter lambda=1 to determine the number of generators
# that are outaged
include("..\\src\\parser.jl")
# network case to model
ps = mp2ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case39.m")# ("..\\data\\case39.m");
# parameters of distributions for line outages and recovery times
#lines_dist = CSV.read("line-distribution-parameters.csv");
s_line = 2.56;#lines_dist[1];
mu_line = 3.66;#lines_dist[2];
sigma_line = 2.43;#lines_dist[3];
# parameters of distributions for generator outages and recovery times
#gens_dist = CSV.read("gen-distribution-parameters.csv");
lambda_gen = 1;#gens_dist[1];
mu_gen = 3.66;#gens_dist[2];
sigma_gen = 2.43;#gens_dist[3];
# number of events to realize the outages and restoration times for
Nevents = 10000;
# number of lines and generators in network case
TotalLines = length(ps.branch[1]);
TotalGens =length(ps.gen[1]);
include("..\\src\\s1-initiate1.jl")
maxLinesOut = TotalLines;
Nlines = zeros(Nevents)
RecovTime = 1;
for n = 1:Nevents
Nlines[n] = init_out_zipf(s_line,maxLinesOut,TotalLines);
#Ngens = init_out_exp(lambda_gen,TotalGens);
RecovTime = [RecovTime; RecoveryTimes(mu_line,sigma_line,Nlines)];
#RecovTimeG = RecoveryTimes(mu_gen,sigma_gen,Ngens);
#lines_state = initiate_state(TotalLines, Nlines);
#gens_state = initiate_state(TotalGens, Ngens);
end
Nline = sort(Nlines);
ccdf_nlines_empir = collect(length(Nline):-1:1)./length(Nline);

RecovTim = RecovTime[2:end];
RecT = sort(RecovTim);
ccdf_rectime_empir = collect(length(RecT):-1:1)./length(RecT);

num = DataFrame(Realizations_of_NumLines = Nlines, Number_Lines_Out = Nline, CCDF = ccdf_nlines_empir);
rec = DataFrame(Realizations_of_RecovTime = RecovTim[:], Recovery_Time = RecT, CCDF = ccdf_rectime_empir);
k = 1:100;
pdf_nlines = (Nevents/zeta(s_line)).*k.^(-s_line);
H_k_s = zeros(k[end]);
for i = 1:k[end]
    for j = 1:i
    H_k_s[i] = H_k_s[i] + 1/(j^s_line);
    end
end
ccdf_lines = (1-H_k_s./zeta(s_line));
NumLinesAnalytic = DataFrame(k=k, pdf_Num_Lines_Out=pdf_nlines, ccdf_Num_Lines_Out=ccdf_lines);
T = 0:0.01:10000;
ccdf_rectime = zeros(length(T));
for i = 1:length(T)
    t = T[i];
    ccdf_rectime[i] = 1-(0.5(1+erf((log.(t)-mu_line)./(sqrt(2)*sigma_line))));
end
pdf_rectime = ((length(RecovTime)-1)*(1/(sigma_line*sqrt(2*pi)))).*T.^(-1).*exp.(-((log.(T)-mu_line).^2).*(1/(2*sigma_line^2)));
RecovTimeAnalytic = DataFrame(T = T, pdf_RecoveryTime=pdf_rectime, ccdf_RecoveryTime=ccdf_rectime)
using Gadfly
NumLines = plot(layer(NumLinesAnalytic, x="k", y="ccdf_Num_Lines_Out", Theme(default_color=colorant"red"), Geom.line), layer(num, x = "Number_Lines_Out", y ="CCDF", Geom.bar), Scale.x_log10);
RestoreT = plot(layer(RecovTimeAnalytic, x="T", y="ccdf_RecoveryTime", Theme(default_color=colorant"red"), Geom.line), layer(rec, x = "Recovery_Time", y = "CCDF", Geom.bar), Scale.x_log10, Scale.y_log10);
draw(PDF("validate-BPA-line-outage.pdf", 10inch, 6inch), NumLines);
draw(PDF("validate-BPA-recovery-time.pdf", 10inch, 6inch), RestoreT);
