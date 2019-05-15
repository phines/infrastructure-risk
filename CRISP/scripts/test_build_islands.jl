#this will be easy
include("..\\src\\CRISP_network.jl")
# pass nodes and link
#make line graph--> test answer is correct


## load the case data
#ps = mp2ps("../data/case6ww.m")
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case6ww\\")

# remove branches
ps.branch[4,:status]=0;
ps.branch[7,:status]=0;
ps.branch[8,:status]=0;
ps.branch[11,:status]=0;

subgraph = find_subgraphs(ps);
ps_islands = build_islands(subgraph,ps)

#test on big system to check for reasonable results
