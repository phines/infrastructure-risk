using CSV; using DataFrames; using SpecialFunctions;
include("CRISP_LSOPF_tests.jl")
include("CRISP_network.jl")

function RLSOPF1!(totalp,ps,failures,recovery_times,Pd_max;t0 = 10, load_cost=0)
    if sum(load_cost)==0
        load_cost = ones(length(ps.shunt.P));
    end
    rec_t = recovery_times[recovery_times.!=0];
    times = sort(rec_t);
    load_shed = zeros(length(times)+2);
    load_shed[1] = 0;
    lines_out = zeros(length(times)+2);
    lines_out[2] = length(failures) - sum(failures);
    # set load shed for the step just before restoration process
    load_shed[2] = sum(load_cost.*(Pd_max - ps.shunt[:P]));
    for i = 1:length(times)
        T = times[i];
        # set failed branches to status 0
        failures[T.>=recovery_times] .= 1;
        lines_out[i+2] = length(failures) - sum(failures);
        # apply to network
        ps.branch[:,:status] = failures;
        #check for islands
        subgraph = find_subgraphs(ps);# add Int64 here hide info here
        M = Int64(findmax(subgraph)[1]);
        ps_islands = build_islands(subgraph,ps);# at some point check for changes in islands and don't run power flows if no change
        ## for every island that changed (eventually)
        for j in 1:M
            psi = ps_subset(ps,ps_islands[j]);
#            if sum(ps_islands[j].bus)==1
#                #Don't need to rerun for single bus islands (already ran any single bus in lsopf).
#                #Note, this will change if generators or loads start getting damaged
#            else
                # run the dcpf
                crisp_dcpf1!(psi);
                # run lsopf
                (dPd, dPg) = crisp_rlopf1(psi,Pd_max[ps_islands[j].shunt]);
                # apply the results
                ps.gen[ps_islands[j].gen,:Pg]  += dPg;
                ps.shunt[ps_islands[j].shunt,:P] += dPd;
                crisp_dcpf1!(psi);
#            end
        end
        # set load shed for this time step
        load_shed[i+2] = sum(load_cost.*(Pd_max - ps.shunt[:P]));
    end
    times = [0.0;t0*1.0;times.+t0*1.0];
    perc_load_served = (sum(load_cost.*Pd_max) .- load_shed)./sum(load_cost.*Pd_max);
    Restore = DataFrame(time = times, load_shed = load_shed, perc_load_served = perc_load_served, num_lines_out = lines_out);
    return Restore
end

function crisp_rlopf1(ps,Pd_max)
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
        sparse(F,F,+Xinv,n,n);
    ### Build the optimization model ###
    m1 = Model(with_optimizer(Clp.Optimizer));
    # variables
    @variable(m1,dPd[1:nd]);
    @variable(m1,dPg[1:ng]);
    @variable(m1,dTheta[1:n]);
    # variable bounds
    @constraint(m1,-Pd.<=dPd.<=(Pd_max./ ps.baseMVA-Pd));
    @constraint(m1,-Pg.<=dPg.<=(ps.gen[:Pmax]-Pg));
    @constraint(m1,dTheta[1] == 0);
    # objective
    @objective(m1,Max,sum(dPd) - 0.9*sum(dPg)); # serve as much load as possible
    # mapping matrix to map loads/gens to buses
    M_D = sparse(D,1:nd,1.0,n,nd);
    M_G = sparse(G,1:ng,1.0,n,ng);
    # Power balance equality constraint
    @constraint(m1,B*dTheta .== M_G*dPg - M_D*dPd);
    # Power flow constraints
    @constraint(m1,-flow_max .<= flow0 + Xinv.*(dTheta[F] - dTheta[T]) .<= flow_max);
    ### solve the model ###
    optimize!(m1);
    # collect/return the outputs
    sol_dPd=value.(dPd);
    sol_dPg=value.(dPg);
    dPd_star = sol_dPd.*ps.baseMVA;
    dPg_star = sol_dPg.*ps.baseMVA;
    return (dPd_star, dPg_star)
end
