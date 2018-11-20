module CRISP
#set up packages
using CSV; using DataFrames; using SpecialFunctions;
include("s1-initiate1.jl")
function line_state(ps,s_line,maxLinesOut,mu_line,sigma_line,orignumLines=0)
# number of lines and generators in network case
TotalLines = length(ps.branch[1]);
Nlines = init_out_zipf(s_line,maxLinesOut,TotalLines);
lines_state = initiate_state(TotalLines, Nlines);
RecovTimeL = RecoveryTimes(mu_line,sigma_line,Nlines);
lines_outage_recovery = RecTime(RecovTimeL,lines_state)
return lines_outage_recovery
end
function gen_state(ps,lambda_gen,mu_gen,sigma_gen,orignumGen)
# number of lines and generators in network case
TotalGens =length(ps.gen[1]);
Ngens = init_out_exp(lambda_gen,TotalGens);
gens_state = initiate_state(TotalGens, Ngens);
RecovTimeG = RecoveryTimes(mu_gen,sigma_gen,Ngens);
gens_outage_recovery = RecTime(RecovTimeG,gens_state)
return gens_outage_recovery
end
end
