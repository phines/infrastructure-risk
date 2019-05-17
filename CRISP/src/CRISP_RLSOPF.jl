using CSV; using DataFrames; using SpecialFunctions; using JuMP; using Clp; using Gurobi
include("CRISP_LSOPF.jl")
include("CRISP_network.jl")

function RLSOPF!(totalp,ps,failures,recovery_times,Pd_max ;t0 = 10, load_cost=0)
    # constants
    tolerance = 1e-6
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
        if M==1
            crisp_dcpf!(ps);
            # run lsopf
            #crisp_rlopf1!(ps,Pd_max,Flow0);
            crisp_rlopf!(ps,Pd_max);
            crisp_dcpf!(ps);
        else
            ps_islands = build_islands(subgraph,ps);# at some point check for changes in islands and don't run power flows if no change
            ## for every island that changed (eventually)

            for j in 1:M
                psi = ps_subset(ps,ps_islands[j]);
                    # run the dcpf
                    crisp_dcpf!(psi);
                    # run lsopf
                    crisp_rlopf!(psi,Pd_max[ps_islands[j].shunt]);
                    # apply the results
                    ps.gen[ps_islands[j].gen,:Pg]  = psi.gen.Pg;
                    ps.shunt[ps_islands[j].shunt,:P] = psi.shunt.P;
                    crisp_dcpf!(psi);
            end
        end
        @assert abs(sum(ps.gen.Pg)-sum(ps.shunt.P))<=tolerance
        # set load shed for this time step
        load_shed[i+2] = sum(load_cost.*(Pd_max - ps.shunt[:P]));
    end
    times = [0.0;t0*1.0;times.+t0*1.0];
    perc_load_served = (sum(load_cost.*Pd_max) .- load_shed)./sum(load_cost.*Pd_max);
    Restore = DataFrame(time = times, load_shed = load_shed, perc_load_served = perc_load_served, num_lines_out = lines_out);
    return Restore
end

function crisp_rlopf!(ps,Pd_max)
    # constants
    tolerance = 1e-6
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    if n>1
        bi = sparse(ps.bus.id,fill(1,n),collect(1:n)) # helps us to find things
        # load data
        nd = size(ps.shunt,1)
        D = bi[ps.shunt[:bus]]
        D_bus = sparse(D,collect(1:nd),1.,n,nd);
        Pd = ps.shunt[:P] ./ ps.baseMVA .* ps.shunt[:status]
        if any(D.<1) || any(D.>n)
            error("Bad indices in shunt matrix")
        end
        #Pd_bus = Array(sparse(D,ones(size(D)),Pd,n,1))
        # gen data
        ng = size(ps.gen,1)
        G = bi[ps.gen[:bus]]
        G_bus = sparse(G,collect(1:ng),1.,n,ng);
        Pg = ps.gen[:Pg] ./ ps.baseMVA .* ps.gen[:status]
        Pg_max = ps.gen[:Pmax] ./ ps.baseMVA .* ps.gen[:status]
        if any(G.<1) || any(G.>n)
            error("Bad indices in gen matrix")
        end
        #Pg_bus = Array(sparse(G,ones(size(G)),Pg,n,1))
        # branch data
        brst = (ps.branch.status.==1)
        F = bi[ps.branch.f[brst]]
        T = bi[ps.branch.t[brst]]
        flow0 = ps.branch.Pf[brst]./ps.baseMVA
        flow_max = ps.branch.rateA[brst]./ps.baseMVA # this could also be rateB
        Xinv = (1 ./ ps.branch.X[brst])
        B = sparse(F,T,-Xinv,n,n) +
            sparse(T,F,-Xinv,n,n) +
            sparse(T,T,+Xinv,n,n) +
            sparse(F,F,+Xinv,n,n)
        ### Build the optimization model ###
        m1 = Model(with_optimizer(Gurobi.Optimizer))
        # variables
        @variable(m1,dPd[1:nd])
        @variable(m1,dPg[1:ng])
        @variable(m1,dTheta[1:n])
        # variable bounds
        @constraint(m1,-Pd.<=dPd.<=((Pd_max./ ps.baseMVA .* ps.shunt[:status]) - Pd));
        @constraint(m1,-Pg.<=dPg.<=Pg_max-Pg);
        @constraint(m1,dTheta[1] == 0);
        # objective
        @objective(m1,Max,sum(dPd)) # serve as much load as possible
        # mapping matrix to map loads/gens to buses
        # Power balance equality constraint
        @constraint(m1,B*dTheta .== G_bus*dPg-D_bus*dPd);#M_G*dPg - M_D*dPd)
        # Power flow constraints
        @constraint(m1,-flow_max .<= flow0 + Xinv.*(dTheta[F] - dTheta[T]) .<= flow_max)
        ### solve the model ###
        optimize!(m1);
        # collect/return the outputs
        sol_dPd=value.(dPd)
        sol_dPg=value.(dPg)
        dPd_star = sol_dPd.*ps.baseMVA
        dPg_star = sol_dPg.*ps.baseMVA
        ps.shunt.P += dPd_star; #changes ps structure
        ps.gen.Pg += dPg_star; #changes ps structure
    else
        if (!isempty(ps.gen) && isempty(ps.shunt)) || (isempty(ps.gen) && !isempty(ps.shunt)) || (isempty(ps.gen) && isempty(ps.shunt))
            ps.gen.Pg  .= ps.gen.Pg.*0.0;
            ps.shunt.P .= ps.shunt.P.*0.0;
        elseif !isempty(ps.gen) && !isempty(ps.shunt)
            Pd = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
            Pg = ps.gen.Pg ./ ps.baseMVA .* ps.gen.status
            Pg_cap = ps.gen.Pmax ./ ps.baseMVA .* ps.gen.status
            if sum(Pg_cap) >= sum(Pd)
                deltaPd = 0.0;
                deltaPg = sum(Pd)-sum(Pg);
            else
                deltaPd = sum(Pg_cap)-sum(Pd);
                deltaPg = sum(Pg_cap)-sum(Pg);
            end
            if length(Pd)  > 1
                deltaPd_star = zeros(length(Pd))
                deltaPd_used = 0;
                for load = 1:length(Pd)
                    if abs(deltaPd) >= abs(deltaPd_used)
                        dload = 0;
                    else
                        if abs(Pd[load]) < abs(deltaPd_used - deltaPd)
                            dload = -Pd[load]
                        else
                            dload = -abs(deltaPd_used - deltaPd)
                        end
                    end
                    deltaPd_star[load] = dload.*ps.baseMVA;
                    deltaPd_used = deltaPd_used + dload;
                end
            else
                deltaPd_star = deltaPd.*ps.baseMVA;
            end
            if length(Pg) > 1
                deltaPg_star = zeros(length(Pg))
                deltaPg_used = 0;
                for gen = 1:length(Pg)
                    deltaPg_star[gen] = deltaPg.*ps.baseMVA;
                    if abs(deltaPg) >= abs(deltaPg_used)
                        dgen = 0;
                    else
                        if abs(Pg[gen]) < abs(deltaPg_used - deltaPg)
                            if deltaPg < 0
                                dgen = -Pg[gen]
                            else
                                dgen = Pg[gen]
                            end
                        elseif deltaPg < 0
                            dgen = -abs(deltaPg_used - deltaPg)
                        else
                            dgen = abs(deltaPg_used - deltaPg)
                        end
                    end
                    deltaPg_star[gen] = dgen.*ps.baseMVA;
                    deltaPg_used = deltaPg_used + dgen;
                end
            else
                deltaPg_star = deltaPg.*ps.baseMVA;
            end
            ps.shunt.P .+= deltaPd_star;
            ps.gen.Pg  .+= deltaPg_star;
        else
            ps.shunt.P = ps.shunt.P.*0.0;
            ps.gen.Pg  = ps.gen.Pg.*0.0;
        end
    end
    #adding criteria that should produce errors if incorrect.
    @assert abs(sum(ps.shunt.P)-sum(ps.gen.Pg))<=2*tolerance
    #@assert 0.<=ps.shunt.P
    return ps
end
function crisp_rlopf1!(ps,Pd_max,Flow0)
    # constants
    tolerance = 1e-6
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
        bi = sparse(ps.bus.id,fill(1,n),collect(1:n)) # helps us to find things
        # load data
        nd = size(ps.shunt,1)
        D = bi[ps.shunt[:bus]]
        D_bus = sparse(D,collect(1:nd),1.,n,nd);
        Pd = ps.shunt[:P] ./ ps.baseMVA .* ps.shunt[:status]
        if any(D.<1) || any(D.>n)
            error("Bad indices in shunt matrix")
        end
        #Pd_bus = Array(sparse(D,ones(size(D)),Pd,n,1))
        # gen data
        ng = size(ps.gen,1)
        G = bi[ps.gen[:bus]]
        G_bus = sparse(G,collect(1:ng),1.,n,ng);
        Pg = ps.gen[:Pg] ./ ps.baseMVA .* ps.gen[:status]
        Pg_max = ps.gen[:Pmax] ./ ps.baseMVA .* ps.gen[:status]
        if any(G.<1) || any(G.>n)
            error("Bad indices in gen matrix")
        end
        #Pg_bus = Array(sparse(G,ones(size(G)),Pg,n,1))
        # branch data
        brst = (ps.branch.status.==1)
        F = bi[ps.branch.f[brst]]
        T = bi[ps.branch.t[brst]]
        flow0 = Flow0[brst]./ps.baseMVA
        flow_max = ps.branch.rateA[brst]./ps.baseMVA # this could also be rateB
        Xinv = (1 ./ ps.branch.X[brst])
        B = sparse(F,T,-Xinv,n,n) +
            sparse(T,F,-Xinv,n,n) +
            sparse(T,T,+Xinv,n,n) +
            sparse(F,F,+Xinv,n,n)
        ### Build the optimization model ###
        m2 = Model(with_optimizer(Clp.Optimizer))
        # variables
        @variable(m2,dPd[1:nd])
        @variable(m2,dPg[1:ng])
        @variable(m2,Theta[1:n])
        # variable bounds
        @constraint(m2,0 .<=dPd.<=(Pd_max./ ps.baseMVA .* ps.shunt[:status]));
        @constraint(m2,0 .<=dPg.<=Pg_max);
        @constraint(m2,Theta[1] == 0);
        # objective
        @objective(m2,Max,sum(dPd)) # serve as much load as possible
        # mapping matrix to map loads/gens to buses
        # Power balance equality constraint
        @constraint(m2,B*Theta .== G_bus*dPg-D_bus*dPd);
        # Power flow constraints
        @constraint(m2,-flow_max .<= Xinv.*(Theta[F] - Theta[T]) .<= flow_max)
        ### solve the model ###
        optimize!(m2);
        # collect/return the outputs
        sol_dPd=value.(dPd)
        sol_dPg=value.(dPg)
        dPd_star = sol_dPd.*ps.baseMVA
        dPg_star = sol_dPg.*ps.baseMVA
        ps.shunt.P += dPd_star; #changes ps structure
        ps.gen.Pg += dPg_star; #changes ps structure

    #adding criteria that should produce errors if incorrect.
    @assert abs(sum(ps.shunt.P)-sum(ps.gen.Pg))<=2*tolerance
    #@assert 0.<=ps.shunt.P
    return ps
end
