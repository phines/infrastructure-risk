#script that tests function find_subgraphs
include("../src/CRISP_network.jl")
include("../src/CRISP_LSOPF_tests.jl")

#test subgraph 1
# load the case data
#ps = mp2ps("../data/case6ww.m")
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case6ww\\")
crisp_dcpf!(ps)
ps0 = deepcopy(ps);
#println(ps)
println(ps.gen.Pg)

#subgraph 1
# remove branches
ps.branch[2,:status]=0;
ps.branch[5,:status]=0;

#find subgraph and print
subgraph = find_subgraphs(ps)
println(subgraph)
#find islands and print
ps_islands = build_islands(subgraph,ps);
ps1 = ps_subset(ps,ps_islands[1]);
crisp_dcpf!(ps1)
crisp_lsopf!(ps1)
println(ps1.gen.Pg)
println(ps1.shunt.P)

#reset network
ps = deepcopy(ps0)
# remove branches
ps.branch[2,:status]=0;
ps.branch[5,:status]=0;
#run load shedding all at once
crisp_dcpf!(ps)
crisp_lsopf!(ps)
println(ps.gen.Pg)
println(ps.shunt.P)
