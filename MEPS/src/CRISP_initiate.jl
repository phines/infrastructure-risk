# uses the distributions of the number of lines outages and the distribution of recovery times
# lines which are fits to BPA data to realize a specifc number of line failures and the recovery
# time of each failure
# it also uses an exponential function with paamter lambda=1 to determine the number of generators
# that are outaged

using SpecialFunctions;
using Random

function line_state!(ps,s_line,maxLinesOut,mu_line,sigma_line;orignumLines=0)
# number of lines and generators in network case
TotalLines = length(ps.branch.f);
Nlines = init_out_zipf_p1(s_line,maxLinesOut,TotalLines);
lines_state = initiate_state(TotalLines, Nlines);
ps.branch.status .= lines_state;
RecovTimeL = RecoveryTimes(mu_line,sigma_line,Nlines);
lines_outage_recovery = RecTime(RecovTimeL,lines_state);
return lines_outage_recovery
end

function line_state2!(ps,s_line,maxLinesOut,mu_line,sigma_line;orignumLines=0)
# number of lines and generators in network case
TotalLines = length(ps.branch[1]);
Nlines = init_out_zipf_a0(s_line,maxLinesOut,TotalLines);
lines_state = initiate_state(TotalLines, Nlines);
ps.branch[:status] = lines_state;
RecovTimeL = RecoveryTimes(mu_line,sigma_line,Nlines);
lines_outage_recovery = RecTime(RecovTimeL,lines_state);
return lines_outage_recovery
end

#No longer using.
# fair sample from distribution that allows 0 outages, but that only picks events with outages
function line_state3!(ps,s_line,maxLinesOut,mu_line,sigma_line;orignumLines=0)
# number of lines and generators in network case
TotalLines = length(ps.branch[1]);
Nlines = init_out_zipf_a0(s_line,maxLinesOut,TotalLines);
while Nlines==0
    Nlines = init_out_zipf_a0(s_line,maxLinesOut,TotalLines);
end
lines_state = initiate_state(TotalLines, Nlines);
ps.branch[:status] = lines_state;
RecovTimeL = RecoveryTimes(mu_line,sigma_line,Nlines);
lines_outage_recovery = RecTime(RecovTimeL,lines_state);
return lines_outage_recovery
end

function gen_state!(ps,lambda_gen,mu_gen,sigma_gen;orignumGen=0)
# number of lines and generators in network case
TotalGens =length(ps.gen.Pg);
Ngens = init_out_exp(lambda_gen,TotalGens);
gens_state = initiate_state(TotalGens, Ngens);
ps.gen.status = gens_state;
ps.gen.Pg = ps.gen.Pg .*gens_state;
RecovTimeG = RecoveryTimes(mu_gen,sigma_gen,Ngens);
gens_outage_recovery = RecTime(RecovTimeG,gens_state)
return gens_outage_recovery
end

function init_out_zipf_p1(s,k,TotalLines;OrigNumLines=TotalLines)
ratioL = TotalLines/OrigNumLines;
# the number of lines outaged probability distribution is fit to a zipf distribution with s = 2.56
# the cdf of a zipf distribution with
H_k_s = zeros(k-1);
for i = 1:k-1
    for j = 1:i
    H_k_s[i] = H_k_s[i] + 1/(j^s);
    end
end
cdf_lines = H_k_s./zeta(s);
P_leqNlinesOut = rand(rng,1);
cdf_lines = H_k_s./zeta(s);
Nlines = sum(P_leqNlinesOut.>=cdf_lines);
Nlines += 1;
Nlines = Int64(round(ratioL*Nlines)) #
return Nlines
end

function init_out_zipf(s,k,TotalLines;OrigNumLines=TotalLines)
ratioL = TotalLines/OrigNumLines;
# the number of lines outaged probability distribution is fit to a zipf distribution with s = 2.56
# the cdf of a zipf distribution with
H_k_s = zeros(k);
for i = 1:k
    for j = 1:i
    H_k_s[i] = H_k_s[i] + 1/(j^s);
    end
end
cdf_lines = H_k_s./zeta(s);
P_leqNlinesOut = rand(rng,1);
cdf_lines = H_k_s./zeta(s);
Nlines = sum(P_leqNlinesOut.>=cdf_lines);
if Nlines==0
    Nlines=1;
end
Nlines = Int64(round(ratioL*Nlines)) #
return Nlines
end

function init_out_zipf_a0(s,k,TotalLines;OrigNumLines=TotalLines)
ratioL = TotalLines/OrigNumLines;
# the number of lines outaged probability distribution is fit to a zipf distribution with s = 2.56
# the cdf of a zipf distribution with
H_k_s = zeros(k);
for i = 1:k
    for n = 1:i
    H_k_s[i] = H_k_s[i] + 1/(n^s);
    end
end
cdf_lines = H_k_s./zeta(s);
P_leqNlinesOut = rand(rng,1);
Nlines = sum(P_leqNlinesOut.>=cdf_lines);
Nlines = Int64(round(ratioL*Nlines)) #
return Nlines
end

function init_out_exp(lambda,TotalGens;OriginalGens=TotalGens)
# since we do not have the outage distribution of generators, I am modeling it temporariliy as
# an exponential random variable with lambda=1. I would expect that the distribution of
# generators that experience outages (not caused by grid dynamics) will be steeper than the
# lines, which is why I am not using the zipf distribution
ratioG = TotalGens/OriginalGens;
Ngens1 = -round.(log.(1 .-rand(rng,1)));
Ngens = Int64(Ngens1[1]);
Ngens = Int64(round(ratioG*Ngens));
if Ngens > TotalGens
    Ngens = TotalGens
end
return Ngens
end

function RecoveryTimes(mu,sigma,N)
# the restoration time of the lines was fit to log normal with underlying_normal_mean=2.66 and
# underlying_normal_variance= 2.43^2
# the cdf of a lognormal with
    N = Int64(N[1]);
    RecovTime = zeros(N);
    cdf_RT_real = rand(rng,N);
    for m = 1:N
        RecovTime[m] = exp(erfinv(2*cdf_RT_real[m]-1)*(sqrt(2)*sigma)).*exp(mu);
    end
    return RecovTime
end

# next two functions take inital number of line or generator outages and picks the lines or
# generators that actually are removed from the network

#THIS IS WHERE WE CAN INTEGRATE DISTRIBUTION OF NUMBER OF HOPS
function initiate_state(TotalS, N);
    Index = collect(1:TotalS);
    State = ones(TotalS);
    Nout = Random.shuffle(rng,Index);
    for i = 1:N
        State[Nout[i]] = 0;
    end
    return State
end

function RecTime(RecovTime,state)
    j = 1;
    RecTime = zeros(length(state));
    for i = 1:length(RecovTime)
        while state[j] == 1
            j = j+1;
        end
        RecTime[j] = RecovTime[i];
        j = j+1;
    end
    data = DataFrame(state = state, recovery_time = RecTime)
    return data
end