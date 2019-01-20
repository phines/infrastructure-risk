#using CRISP_LSOPF
using CSV
include("../src/CRISP_RLSOPF")
include("../src/CRISP_LSOPF1.jl")
include("../src/parser.jl")

# load the case data
ps = mp2ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case39.m") #("../data/case6ww.m")
#ps = mp2ps("../data/case6ww.m")
crisp_dcpf!(ps)
ps0 = deepcopy(ps)
# remove branches
s_line = 2.56;#lines_dist[1];
maxLinesOut = 70;
mu_line = 3.66;#lines_dist[2];
sigma_line = 2.43;#lines_dist[3];
# parameters of distributions for generator outages and recovery times
#gens_dist = CSV.read("gen-distribution-parameters.csv");
lambda_gen = 1;#gens_dist[1];
mu_gen = 3.66;#gens_dist[2];
sigma_gen = 2.43;#gens_dist[3];

#add to ends of function in CRISP-electricity.jl when original number of
# lines and gens are known.
orignumLines = 0;
orignumGen = 0;

include("..\\src\\CRISP-electricity.jl")
# step 1
Lines_Init_State = line_state(ps,s_line,maxLinesOut,mu_line,sigma_line,orignumLines)
CSV.write("step1_out.csv",Lines_Init_State)
print(Lines_Init_State)
lines_status=Lines_Init_State[1];
#lines_status = ones(length(lines_stat));#zeros(length(lines_stat));
#lines_status[1]=0;
#lines_status[5]=0;
for i=1:length(lines_status)
    if lines_status[i] == 0
        ps.branch[i,:status]=0;
    end
end
# run the dcpf
crisp_dcpf!(ps)
# run lsopf
(dPd, dPg) = crisp_lsopf(ps)
# apply the results
ps.gen[:Pg]  += dPg
ps.shunt[:P] += dPd
crisp_dcpf!(ps)
load_shed[i] = sum(ps.shunt[:P]);
