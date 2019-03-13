#using CRISP_LSOPF
using CSV
include("../src/CRISP_initiate.jl")
include("../src/CRISP_LSOPF.jl")
include("../src/CRISP_network.jl")

# load the case data
#ps = import_ps("../data/case39/")
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case39\\")
#ps = mp2ps("../data/case6ww.m")
# run the dcpf
crisp_dcpf!(ps)
ps0 = deepcopy(ps)

#define parameters for removing branches
s_line = 2.56;#lines_dist[1];
maxLinesOut = 70;
mu_line = 3.66;#lines_dist[2];
sigma_line = 2.43;#lines_dist[3];

# step 1
Lines_Init_State = line_state!(ps,s_line,maxLinesOut,mu_line,sigma_line);
CSV.write("results\\test_step1_out.csv",Lines_Init_State)
print(Lines_Init_State)
failures=Lines_Init_State[1];
ps.branch.status = failures;

#check for islands
subgraph = find_subgraphs(ps);
M = findmax(subgraph)[1];
ps_islands = build_islands(subgraph,ps);
for i in 1:M
    psi = ps_subset(ps,ps_islands[i]);
    ## run step 2
    crisp_lsopf!(psi)
    add_changes!(ps,psi,ps_islands[i]);
    crisp_dcpf!(psi)
end

println("Before");
print(sum(ps0.gen[:Pg]));
println(" ")
println("After");
print(sum(ps.gen[:Pg]));
