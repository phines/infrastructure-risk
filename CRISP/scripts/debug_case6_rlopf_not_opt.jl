using CSV
include("..\\src\\CRISP_initiate.jl")
include("..\\src\\CRISP_LSOPF.jl")
include("..\\src\\CRISP_RLSOPF.jl")
include("..\\src\\CRISP_RT.jl")
include("..\\src\\CRISP_network.jl")


ps = import_ps("data\\case6ww\\")
crisp_dcpf!(ps)
total = sum(ps.shunt[:P]);
Pd_max = deepcopy(ps.shunt[:P]);
ps0 = deepcopy(ps);
Flow0 = ps.branch.Pf;
Lines_Init_State = CSV.read("C:\\Users\\mkellygo\\Documents\\GitHub\\infrastructure-risk\\CRISP\\results\\experiments\\7\\res_out_case6ww_p1 IC916.csv")
state = Lines_Init_State[:,1];
recovery_times = Lines_Init_State[:,2];
failures = state;
ps.branch.status = failures;
    #check for islands
    subgraph = find_subgraphs(ps);
    M = Int64(findmax(subgraph)[1]);
    if M>1
        ps_islands = build_islands(subgraph,ps);
        for i in 1:M
            psi = ps_subset(ps,ps_islands[i]);
            # run the dcpf
            crisp_dcpf!(psi);
            # run lsopf
            crisp_lsopf!(psi);
            ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
            ps.shunt[ps_islands[i].shunt,:P] = psi.shunt.P
            crisp_dcpf!(psi);
        end
            @assert 10^(-6)>=abs(sum(ps.shunt.P)-sum(ps.gen.Pg))
            @assert total>=sum(ps.shunt.P)
    else
        crisp_dcpf!(ps);
        crisp_lsopf!(ps);
        crisp_dcpf!(ps);
    end
    @assert total>=sum(ps.shunt.P)
    @assert 10^(-6)>=abs(sum(ps.shunt.P)-sum(ps.gen.Pg))
    tolerance = 10^(-10);
    ## run step 3
    Restore = RLSOPF!(total,ps,failures,recovery_times,Pd_max);#,load_cost) # data frame [times, load shed in cost per hour]
    ## run step 4
    ResilienceTri = crisp_res(Restore);
