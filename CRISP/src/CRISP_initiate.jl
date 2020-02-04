# uses the distributions of the number of lines outages and the distribution of recovery times
# lines which are fits to BPA data to realize a specifc number of line failures and the recovery
# time of each failure
# it also uses an exponential function with paamter lambda=1 to determine the number of generators
# that are outaged

using SpecialFunctions
using Random
using CSV
using DataFrames
include("CRISP_network.jl")
include("CRISP_LSOPF.jl")


function Outages(Num,ps_folder;param_file = "",cascade=true,comms=true)
    #constants
    debug=1;
    tolerance1 = 10^(-6);
    mult_factor = 1.5
    ## Num = number of failure scenarios to run through
    # initialize vector of costs from events
    ## load the case data
    ps = import_ps("$ps_folder")
    ps.shunt = ps.shunt[ps.shunt.P .!=0.0,:]
    crisp_dcpf_g1_s!(ps)
    if isempty(param_file)
        # parameters of distributions for line outages and recovery times
        s_line = 2.56;#lines_dist[1];
        maxLinesOut = length(ps.branch.f); # => k in zipf distribution
        mu_line = 3.66;#lines_dist[2];
        sigma_line = 2.43;#lines_dist[3];
        lambda_gen = 1;
    end
    for iterat in 1:Num
        # step 1
        if cascade
            Lines_Init_State = line_state_cascade!(ps,s_line,maxLinesOut,mu_line,sigma_line)
            if comms
                dt = 10
                subgraph = find_subgraphs(ps);
                M = Int64(findmax(subgraph)[1]);
                ps_islands = build_islands(subgraph,ps);
                for i in 1:M
                    psi = ps_subset(ps,ps_islands[i]);
                    # run lsopf
                    crisp_lsopf_g1_s!(psi,dt);
                    ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
                    ps.shunt[ps_islands[i].shunt,:status] = psi.shunt.status
                    ps.storage[ps_islands[i].storage,:E] = psi.storage.E
                    ps.storage[ps_islands[i].storage,:Ps] = psi.storage.Ps
                end
                line_rec_times_add_comms!(ps,Lines_Init_State,mult_factor)
            end
        else
            Lines_Init_State = line_state!(ps,s_line,maxLinesOut,mu_line,sigma_line)
        end
        Gens_Init_State = gen_state!(ps,lambda_gen,mu_line,sigma_line)
        if debug==1
            CSV.write("data/outage_data/communication_factor/out_case73_noPWS_lx2_n-1_lines$iterat.csv", Lines_Init_State)
            CSV.write("data/outage_data/communication_factor/out_case73_noPWS_lx2_n-1_gens$iterat.csv", Gens_Init_State)
        end
    end
end

function Outages_ss(Num,ps_folder,out_folder,outfile,bins_lines,bins_gens,gtrip;cascade=true,param_file = "")
    debug=1;
    tolerance1 = 10^(-6);
    nlines = length(bins_lines)
    ngens = length(bins_gens)
    ## load the case data
    ps = import_ps("$ps_folder")
    ps.shunt = ps.shunt[ps.shunt.P .!=0.0,:]
    if isempty(param_file)
        # parameters of distributions for line outages and recovery times
        s_line = 2.56;#lines_dist[1];
        maxLinesOut = length(ps.branch.f); # => k in zipf distribution
        mu_line = 3.66;#lines_dist[2];
        sigma_line = 2.43;#lines_dist[3];
        lambda_gen = 1;
    end
    TotalLines = length(ps.branch.f);
    diameter = find_diameter(ps);
    for iterat in 1:Num
        for N in 1:nlines
            
            # step 1
            if cascade
                lines_state = cascade!(ps, TotalLines, N, diameter);
                RecovTimeL = RecoveryTimes(mu_line, sigma_line, N);
                Lines_Init_State = RecTime(RecovTimeL, lines_state);
            else
                lines_state = initiate_state(TotalLines, N);
                RecovTimeL = RecoveryTimes(mu_line,sigma_line,N);
                Lines_Init_State = RecTime(RecovTimeL,lines_state);
            end
            if gtrip
                line_state = Lines_Init_State.state
                Gens_Init_State = gen_trip!(ps,line_state,mu_line,sigma_line)
                if debug==1
                    if isdir(out_folder*"/$N")
                    else mkdir(out_folder*"/$N") end
                    CSV.write(out_folder*"/$N/"*outfile*"_lines$iterat.csv", Lines_Init_State)
                    CSV.write(out_folder*"/$N/"*outfile*"_gens$iterat.csv", Gens_Init_State)
                end
            else
                for G in 1:ngens
                    TotalGens = length(ps.gen.bus);
                    gens_state = initiate_state(TotalGens,G);
                    RecovTimeG = RecoveryTimes(mu_line,sigma_line,G);
                    Gens_Init_State = RecTime(RecovTimeG,gens_state)
                    if debug==1
                        if isdir(out_folder*"/$N-$G")
                        else mkdir(out_folder*"/$N-$G") end
                        CSV.write(out_folder*"/$N-$G/"*outfile*"_lines$iterat.csv", Lines_Init_State)
                        CSV.write(out_folder*"/$N-$G/"*outfile*"_gens$iterat.csv", Gens_Init_State)
                    end
                end
            end
        end
    end
end

function make_bins_pb(bin_l, bin_h,zipf,geo;s=2.56,lambda=1)
    # the number of lines outaged probability distribution is fit to a zipf distribution with s = 2.56
    # the cdf of a zipf distribution
    if zipf = true
        k = length(bin_l:bin_h)
        H_k_s = zeros(k-1);
        for i = bin_l:(bin_h-1)
            for j = 1:i
            H_k_s[i] = H_k_s[i] + 1/(j^s);
            end
        end
        cdf_lines = H_k_s./zeta(s);
        s = 1/sum(cdf_lines);
        scaled_cdf = cdf_lines.*s;
        P_leqNlinesOut = rand(rng,1);
        Nlines = sum(P_leqNlinesOut.>=cdf_lines);
        Nlines += 1;
        Nlines = Int64(round(ratioL*Nlines)) #
        println(Nlines)
        return Nlines
    elseif geo = true
    else
    end
end

function line_state_cascade!(ps,s_line,maxLinesOut,mu_line,sigma_line;orignumLines=0)
# number of lines and generators in network case
TotalLines = length(ps.branch.f);
Nlines = init_out_zipf(s_line, maxLinesOut, TotalLines);
diameter = find_diameter(ps)
lines_state = cascade!(ps, TotalLines, Nlines, diameter);
RecovTimeL = RecoveryTimes(mu_line, sigma_line, Nlines);
lines_outage_recovery = RecTime(RecovTimeL, lines_state);
return lines_outage_recovery
end

function line_state!(ps,s_line,maxLinesOut,mu_line,sigma_line;orignumLines=0)
# number of lines and generators in network case
TotalLines = length(ps.branch.f);
Nlines = init_out_zipf(s_line,maxLinesOut,TotalLines);
lines_state = initiate_state(TotalLines, Nlines);
ps.branch.status .= lines_state;
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

function gen_trip!(ps,line_state,mu_gen,sigma_gen;orignumGen=0)
# find generators that trip due to line outages in network case
gens_state = initiate_trips(ps,line_state);
Ngens = length(gens_state) - sum(gens_state);
ps.gen.status = gens_state;
ps.gen.Pg = ps.gen.Pg .*gens_state;
RecovTimeG = RecoveryTimes(mu_gen,sigma_gen,Ngens);
gens_outage_recovery = RecTime(RecovTimeG,gens_state)
return gens_outage_recovery
end

function initiate_trips(ps,line_state;hops=1);
# useful things
tolerance = 10^(-6);
Ng = length(ps.gen.bus);
Nl = length(ps.branch.f);
Pg = ps.gen.Pg;
F = ps.branch.f;
T = ps.branch.t;
#set output
gen_state = ones(Ng);
#for i in 1:hops
    for g in 1:Ng
        if abs(Pg[g]) > tolerance
            b = ps.gen.bus[g];
            lf = F .== b;
            lt = T .== b;
            Total = sum(lf) + sum(lt);
            Num = (length(line_state[lf])+length(line_state[lt]) - (sum(line_state[lf])+sum(line_state[lt])))
            if rand(rng) <= Num/Total
                gen_state[g] = 0
            end
        end
    end
#end
return gen_state
end

function init_out_zipf(s,k,TotalLines;OrigNumLines=TotalLines)
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
println(Nlines)
return Nlines
end

function init_out_exp(lambda,TotalGens;OriginalGens=TotalGens)
# since we do not have the outage distribution of generators, I am modeling it temporariliy as
# an exponential random variable with lambda=1. I would expect that the distribution of
# generators that experience outages (not caused by grid dynamics) will be steeper than the
# lines, which is why I am not using the zipf distribution
ratioG = TotalGens/OriginalGens;
Ngens1 = -floor.(log.(1 .-rand(rng,1)));
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

function cascade!(ps, TotalLines, Nlines, diameter);
    #hard coded pdf of cascade line distance in cascade
    CascadeHopsData = CSV.File("data/cascade_data/LineDistanceFreq.csv")  |> DataFrame
    cdf = CascadeHopsData.cdf
    distance = CascadeHopsData.distance
    # randomly pick initial line outage
    init_line = rand(rng,1:TotalLines)
    # define output vector of line states
    lines_status = zeros(TotalLines)
    if Nlines == 1
        lines_status .= 1
        lines_status[init_line] = 0
    elseif Nlines == TotalLines
    elseif Nlines == (TotalLines-1)
        lines_status[init_line] = 1
    else
        lines_status .= 1
        lines_status[init_line] = 0
        for n in 1:(Nlines-1)
            matchPr = (rand(rng,1).>=cdf)
            if sum(matchPr) == 0
                hop = 1
            else
                hop = distance[matchPr][end]
            end
            if hop > diameter
                hop = diameter
            end
            new_line = [];
            c = 0
            while isempty(new_line) && c<=10
                c = c+1
                new_line = find_lines_n_hops(ps,lines_status,hop)
                if isempty(new_line)
                    new_Pr = rand(rng,1)
                    matchPr = new_Pr .>= cdf
                    if sum(matchPr) .== 0
                        hop = 1
                    else
                        hop = distance[matchPr][end]
                    end
                    if hop > diameter
                        hop = diameter
                    end
                end
            end
            if isempty(new_line)
                Index = 1:length(lines_status)
                new_line = Index[lines_status .== 1][1]
            end
            lines_status[new_line] = 0
        end
    end
    ps.branch.status[lines_status.==0] .= 0;
    return lines_status
end

function line_rec_times_add_comms!(ps,Lines_Init_State,mult_factor)
    line_status  = Lines_Init_State[:,1]
    r_time = Lines_Init_State[:,2]
    nl = size(ps.branch.f)[1]
    r_time_new = zeros(nl)
    for b in 1:nl
        if line_status[b] .!= 1
            fr = ps.branch.f[b]
            to = ps.branch.t[b]
            if (sum(fr .== ps.shunt.bus) != 0) || (sum(to .== ps.shunt.bus) != 0)
                if sum(ps.shunt.status[ps.shunt.bus .== fr] .!= 1) != 0
                    r_time_new[b] = r_time[b].*mult_factor
                    println("comms effect recovery time of line $b")
                elseif sum(ps.shunt.status[ps.shunt.bus .== to] .!= 1) != 0
                    r_time_new[b] = r_time[b].*mult_factor
                    println("comms effect recovery time of line $b")
                else
                    r_time_new[b] = r_time[b]
                end
            else
                r_time_new[b] = r_time[b]
            end
        end
    end
    Lines_Init_State = DataFrame(state = line_status, recovery_time = r_time_new)
    return ps, Lines_Init_State
end

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
