#using CRISP_LSOPF
using CSV
#include code for all necessary steps (2,3,4,and grid segmenting)
include("..\\src\\CRISP_LSOPF_1.jl")
include("..\\src\\CRISP_RLSOPF.jl")
include("..\\src\\CRISP_RT.jl")
include("..\\src\\CRISP_network_segments.jl")
include("..\\src\\parser.jl")
#include("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\src\\parser.jl")
## load the case data
ps = mp2ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case6ww.m") #case39.m")
#ps = mp2ps("../data/case6ww.m")
crisp_dcpf!(ps)
total = sum(ps.shunt[:P]);
Pd_max = deepcopy(ps.shunt[:P]);

# remove branches that result in islands
ps.branch[4,:status]=0;
ps.branch[8,:status]=0;

#check for islands
subgraph = find_subgraphs(ps);
M = Int64(findmax(subgraph)[1])
if  M >=2
    ps_islands = build_islands(subgraph,ps)
    for i in 1:M
        psi = ps_subset(ps,ps_islands[i])
        ## run step 2
        # run the dcpf
        crisp_dcpf_islands!(ps,ps_islands)
        # run lsopf
        (dPd, dPg) = crisp_lsopf_ilands(ps,ps_islands)
        # apply the results
        ps.gen[:Pg]  += dPg
        ps.shunt[:P] += dPd
        crisp_dcpf_islands!(ps,ps_islands)
    end
    ## run step 3
    Restore = RLSOPF!(total,ps,failures,recovery_times,Pd_max)#,load_cost) # data frame [times, load shed in cost per hour]

    ## run step 4
    ResilienceTri = crisp_res(Restore);
else
    ## run step 2
    # run the dcpf
    crisp_dcpf!(ps)
    # run lsopf
    (dPd, dPg) = crisp_lsopf(ps)
    # apply the results
    ps.gen[:Pg]  += dPg
    ps.shunt[:P] += dPd
    crisp_dcpf!(ps)

    ## run step 3
    Restore = RLSOPF!(total,ps,failures,recovery_times,Pd_max)#,load_cost) # data frame [times, load shed in cost per hour]

    ## run step 4
    ResilienceTri = crisp_res(Restore);
end
