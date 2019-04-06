#script that tests function find_subgraphs
include("../src/CRISP_network.jl")
include("../src/CRISP_LSOPF_tests.jl")

#test subgraph 1
# load the case data
#ps = mp2ps("../data/case6ww.m")
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case6ww\\")
crisp_dcpf1!(ps)
ps0 = deepcopy(ps);
#println(ps)
println(ps.gen.Pg)

#subgraph 1
# remove branches
ps.branch[1,:status]=0;
ps.branch[2,:status]=0;
ps.branch[3,:status]=0;
ps.branch[4,:status]=0;
ps.branch[5,:status]=0;
ps.branch[7,:status]=0;
ps.branch[8,:status]=0;
ps.branch[11,:status]=0;
#find subgraph and print
subgraph = find_subgraphs(ps)
println(subgraph)
#find islands and print
ps_islands = build_islands(subgraph,ps);

println(ps_islands[2])
println(ps_islands[3])
ps1 = ps_subset(ps,ps_islands[1]);
println(ps_islands[1])
println(ps1)
crisp_dcpf1!(ps1)
crisp_lsopf1!(ps1)

println(ps1.gen.Pg)
println(ps1.shunt.P)
ps2 = ps_subset(ps,ps_islands[2]);
crisp_dcpf1!(ps2)
crisp_lsopf1!(ps2)
println(ps2.gen.Pg)
println(ps2.shunt.P)
ps3 = ps_subset(ps,ps_islands[3]);
crisp_dcpf1!(ps3)
crisp_lsopf1!(ps3)
println(ps3.gen.Pg)
println(ps3.shunt.P)

crisp_lsopf1!(ps)
println(ps.gen.Pg)
