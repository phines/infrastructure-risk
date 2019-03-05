#set up packages
using CSV; using DataFrames; using SpecialFunctions;

include("CRISP_LSOPF_1.jl")
include("CRISP_network_segments.jl")

function RLSOPF!(totalp,ps,failures,recovery_times,Pd_max;load_cost=0)
    if load_cost==0
        load_cost = ones(length(ps.shunt[:P]));
    end
    rec_t = recovery_times[recovery_times.!=0];
    times = sort(rec_t);
    load_shed = zeros(length(times)+1);
    # set load shed for the step just before restoration process
    load_shed[1] = sum(load_cost.*(Pd_max - ps.shunt[:P]));
    for i = 1:length(times)
        T = times[i];
        # set failed branches to status 0
        failures[T.>=recovery_times] .= 1;
        # apply to network
        ps.branch[:,:status] = failures;
        #check for islands
        subgraph = find_subgraphs(ps);
        M = Int64(findmax(subgraph)[1]);
        ps_islands = build_islands(subgraph,ps);
        for i in 1:M
            psi = ps_subset(ps,ps_islands[i]);
            # run the dcpf
            crisp_dcpf!(psi);
            # run lsopf
            (dPd, dPg) = crisp_rlopf(psi,Pd_max[ps_islands[i].shunt]);
            # apply the results
            ps.gen[ps_islands[i].gen,:Pg]  += dPg;
            ps.shunt[ps_islands[i].shunt,:P] += dPd;
            crisp_dcpf!(psi);
        end
        # set load shed for this time step
        load_shed[i+1] = sum(load_cost.*(Pd_max - ps.shunt[:P]));
    end
    times = [0;times];
    perc_load_served = (sum(load_cost.*Pd_max) - load_shed)./sum(load_cost.*Pd_max);
    Restore = DataFrame(time = times, load_shed = load_shed, perc_load_served = perc_load_served);
    return Restore
end

function crisp_rlopf(ps,Pd_max)
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1); # the number of buses
    bi = sparse(ps.bus[:id],fill(1,n),collect(1:n)); # helps us to find things
    # load data
    nd = size(ps.shunt,1);
    D = bi[ps.shunt[:bus]];
    Pd = ps.shunt[:P] ./ ps.baseMVA .* ps.shunt[:status];
    # gen data
    ng = size(ps.gen,1);
    G = bi[ps.gen[:bus]];
    Pg = ps.gen[:Pg] ./ ps.baseMVA .* ps.gen[:status];
    # branch data
    brst = (ps.branch[:status].==1);
    F = bi[ps.branch[brst,:f]];
    T = bi[ps.branch[brst,:t]];
    flow0 = ps.branch[brst,:Pf]./ps.baseMVA;
    flow_max = ps.branch[brst,:rateA]./ps.baseMVA; # this could also be rateB
    Xinv = (1 ./ ps.branch[brst,:X])
    B = sparse(F,T,-Xinv,n,n) +
        sparse(T,F,-Xinv,n,n) +
        sparse(T,T,+Xinv,n,n) +
        sparse(F,F,+Xinv,n,n)
    ### Build the optimization model ###
    m1 = Model(solver = ClpSolver())
    # variables
    @variable(m1,dPd[1:nd])
    @variable(m1,dPg[1:ng])
    @variable(m1,dTheta[1:n])
    # variable bounds
    @constraint(m1,-Pd.<=dPd.<=(Pd_max./ ps.baseMVA-Pd))
    @constraint(m1,-Pg.<=dPg.<=(ps.gen[:Pmax]-Pg))
    @constraint(m1,dTheta[1] == 0)
    # objective
    @objective(m1,Max,sum(dPd) - 0.9*sum(dPg)) # serve as much load as possible
    # mapping matrix to map loads/gens to buses
    M_D = sparse(D,1:nd,1.0,n,nd);
    M_G = sparse(G,1:ng,1.0,n,ng);
    # Power balance equality constraint
    @constraint(m1,B*dTheta .== M_G*dPg - M_D*dPd)
    # Power flow constraints
    @constraint(m1,-flow_max .<= flow0 + Xinv.*(dTheta[F] - dTheta[T]) .<= flow_max)
    ### solve the model ###
    solve(m1)
    # collect/return the outputs
    dPd_star = getvalue(dPd).*ps.baseMVA;
    dPg_star = getvalue(dPg).*ps.baseMVA;
    return (dPd_star, dPg_star)
end
