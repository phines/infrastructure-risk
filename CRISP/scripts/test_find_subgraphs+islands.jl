#script that tests function find_subgraphs
include("../src/CRISP_network.jl")

#test subgraph 1
# load the case data
#ps = mp2ps("../data/case6ww.m")
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case6ww\\")
ps0 = deepcopy(ps);
# remove branches
ps.branch[4,:status]=0;
ps.branch[7,:status]=0;
ps.branch[8,:status]=0;
ps.branch[11,:status]=0;
#find subgraph and print
subgraph = find_subgraphs(ps)
println(subgraph)
#find islands and print
ps_islands = build_islands(subgraph,ps);
println(ps_islands[1])
println(ps_islands[2])

#subgraph 2
# load the case data
ps = deepcopy(ps0);
# remove branches
ps.branch[4,:status]=0;
ps.branch[7,:status]=0;
#find subgraph and print
subgraph = find_subgraphs(ps);
println(subgraph)
#find islands and print
ps_islands = build_islands(subgraph,ps);
println(ps_islands[1])

#subgraph 3
# load the case data
ps = deepcopy(ps0);
# remove branches
ps.branch[11,:status]=0;
#find subgraph and print
subgraph = find_subgraphs(ps)
println(subgraph)
#find islands and print
ps_islands = build_islands(subgraph,ps);
println(ps_islands[1])

#subgraph 4
# load the case data
ps = deepcopy(ps0);
# remove branches
ps.branch[1,:status]=0;
#find subgraph and print
subgraph = find_subgraphs(ps)
println(subgraph)
#find islands and print
ps_islands = build_islands(subgraph,ps);
println(ps_islands[1])
