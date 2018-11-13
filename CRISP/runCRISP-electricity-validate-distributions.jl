#set up packages
using CSV; using DataFrames; using SpecialFunctions;

# uses the distributions of the number of lines outages and the distribution of recovery times
# lines which are fits to BPA data to realize a specifc number of line failures and the recovery
# time of each failure
# it also uses an exponential function with paamter lambda=1 to determine the number of generators
# that are outaged
include("parser.jl")
# network case to model
ps = mp2ps("case39.m");
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
include("s1-initiate1.jl")
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
RecovTim = RecovTime[2:end];
RecTime = sort(RecovTim);
cdf_rectime_empir = (length(RecTime):-1:1)./length(RecTime)

num = DataFrame(Realizations_of_NumLines = Nlines);
rec = DataFrame(Realizations_of_RecovTime = RecovTim[:]);
k = 1:50;
pdf_nlines = (Nevents/zeta(s_line)).*k.^(-s_line);
NumLinesAnalytic = DataFrame(k=k, pdf_Num_Lines_Out=pdf_nlines);
t = 0:0.001:500
cdf_rectime = 0.5(1+erf((log.(t)-mu_line)./(sqrt(2)*sigma_line)));
pdf_rectime = ((length(RecovTime)-1)*(1/(sigma_line*sqrt(2*pi)))).*t.^(-1).*exp.(-((log.(t)-mu_line).^2).*(1/(2*sigma_line^2)));
RecovTimeAnalytic = DataFrame(t = t, pdf_RecoveryTime=pdf_rectime, cdf_RecoveryTime=cdf_rectime)
using Gadfly
var_LoadsFulA = plot(layer(NumLinesAnalytic, x="k", y="pdf_Num_Lines_Out", Theme(default_color=colorant"red"), Geom.line), layer(num, x = "Realizations_of_NumLines", Geom.histogram), Scale.x_log10, Scale.y_log10);
var_LoadsFulV = plot(layer(RecovTimeAnalytic, x="t", y="pdf_RecoveryTime", Theme(default_color=colorant"red"), Geom.line), layer(rec, x = "Realizations_of_RecovTime", Geom.histogram), Scale.x_log10);
draw(PDF("validate-BPA-distributions2.pdf", 10inch, 6inch), hstack(var_LoadsFulA, var_LoadsFulV));
