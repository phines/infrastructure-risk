include("..\\src\\CRISP_initiate.jl")
include("..\\src\\CRISP_LSOPF_gen.jl")
include("..\\src\\CRISP_RLSOPF_gen.jl")
include("..\\src\\CRISP_RT.jl")
include("..\\src\\CRISP_network_gen.jl")

function secure_ps!(ps)
    ps0 = deepcopy(ps);
    n = length(ps.branch.Pf);
    maxflow = zeros(n);
    flows = zeros(n);
    for j = 1:n
        println(j)
        ps = deepcopy(ps0);
        ps.branch.status[j] = 0;
        #subgraph = find_subgraphs(ps);# add Int64 here hide info here
        #M = Int64(findmax(subgraph)[1]);
        #if M==1
            ps1 = crisp_dcpf_g!(ps);
            crisp_dcpf_g!(ps)
            flows = abs.(ps.branch.Pf);
            maxflow = (maximum([flows maxflow], dims=2));
        #else
            #ps_islands = build_islands(subgraph,ps);# at some point check for changes in islands and don't run power flows if no change
            #for j in 1:M
                #psi = ps_subset(ps,ps_islands[j]);
                    # run the dcpf
                    #crisp_dcpf!(psi);
                    #flows[ps_islands[j].branch] .= abs.(psi.branch.Pf)
            #end
            #maxflow = (maximum([flows maxflow], dims=2));
        #end
    end
    overloads = maxflow .> ps.branch.rateA
    ps.branch.rateA[overloads[:,1]] .= round.(maxflow[overloads[:,1]].*1.5);
    return ps
end
