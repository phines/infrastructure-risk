using CSV;
using DataFrames;
using SpecialFunctions;
using JuMP;
using Clp;
using Gurobi;
using Cbc;
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
        load_shed[i+2] = sum(load_cost.*(Pd_max - ps.shunt.P));
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
        bi = sparse(ps.bus.id,fill(1,n),collect(1:n));
        # load data
        nd = size(ps.shunt,1)
        D = bi[ps.shunt.bus]
        D_bus = sparse(D,collect(1:nd),1.,n,nd);
        Pd1 = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
        if any(D.<1) || any(D.>n)
            error("Bad indices in shunt matrix")
        end
        # gen data
        ng = size(ps.gen.Pg,1)
        G = bi[ps.gen.bus]
        G_bus = sparse(G,collect(1:ng),1.,n,ng);
        Pg1 = ps.gen.Pg ./ ps.baseMVA .* ps.gen.status
        ug1 = ones(ng); ug1[Pg1.==0] .= 0;
        #Pg = ps.gen[gst,:Pg] ./ ps.baseMVA .* ps.gen[gst,:status]
        Pg_max = ps.gen.Pmax ./ ps.baseMVA .* ps.gen.status
        Pg_min = ps.gen.Pmin ./ ps.baseMVA .* ps.gen.status
        if any(G.<1) || any(G.>n)
            error("Bad indices in gen matrix")
        end
        m = Model(with_optimizer(Gurobi.Optimizer))
        # variables
        @variable(m, Pd[1:nd]) # demand
        @variable(m, Pg[1:ng]) # generation
        # variable bounds constraints
        @constraint(m, 0.0 .<= Pd .<= Pd_max./ ps.baseMVA .* ps.shunt[:status]) # load served limits
        @constraint(m, 0 .<= Pg .<= Pg_max) # generator power limits
        # power balance
        @constraint(m, 0 .== G_bus*Pg-D_bus*Pd)
        # objective
        @objective(m, Max, 100*sum(sum(Pd)));
        ## SOLVE! ##
        optimize!(m)
        sol_Pd=value.(Pd)
        sol_Pg=value.(Pg)
        dPd_star = sol_Pd.*ps.baseMVA
        dPg_star = sol_Pg.*ps.baseMVA
        ps.shunt.P .== dPd_star; #changes ps structure
        ps.gen.Pg .== dPg_star; #changes ps structure
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

## CRISP_RLSOPF_gen_orig.jl
function RLSOPF_g_orig!(ps,l_failures,g_failures,l_recovery_times,g_recovery_times,gen_startup,Pd_max;t0 = 10, load_cost=0)
    #add columns to keep track of the time each generator is on or off
    if sum(names(ps.gen).==:time_on) == 0
        ps.gen.time_on = zeros(length(ps.gen.Pg));
    end
    if sum(names(ps.gen).==:time_off) == 0
        ps.gen.time_off = zeros(length(ps.gen.Pg));
    end
    # constants
    deltaT = 5; # time step in minutes;
    tolerance = 1e-6
    if sum(load_cost)==0
        load_cost = ones(length(ps.shunt.P));
    end
    l_rec_t = l_recovery_times[l_recovery_times.!=0];
    g_rec_t = g_recovery_times[g_failures.==0] + gen_startup[g_failures.==0];
    g_rec = zeros(length(g_failures));
    g_rec[g_failures.==0] = g_rec_t;
    Time = sort([l_rec_t; g_rec_t])
    total_i = length(Time);
    load_shed = zeros(total_i+2);
    load_shed[1] = 0;
    lines_out = zeros(total_i+2);#zeros(length(l_times)+2);
    lines_out[2] = length(l_failures) - sum(l_failures);
    gens_out = zeros(total_i+2);#zeros(length(l_times)+2);
    gens_out[2] = length(g_failures) - sum(g_failures);
    # set load shed for the step just before restoration process
    load_shed[2] = sum(load_cost.*(Pd_max - ps.shunt[:P]));
    time_gen_off = ones(length(ps.gen.Pg));
    time_gen_off[ps.gen.Pg .!=0] .=0;
    for i = 1:length(Time)
        T = Time[i];
        # set failed branches to status 0
        l_failures[T.>=l_recovery_times] .= 1;
        lines_out[i+2] = length(l_failures) - sum(l_failures);
        # apply to network
        ps.branch[:,:status] = l_failures;
        # set newly available generators to status 1;
        current_gen_out = ps.gen.status;
        # update generators operational
        g_failures[T.>=g_rec] .= 1;
        ps.gen.status = g_failures;
        gens_out[i+2] = length(g_failures) - sum(g_failures);
        #check for islands
        subgraph = find_subgraphs(ps);# add Int64 here hide info here
        M = Int64(findmax(subgraph)[1]);
        if M==1
            crisp_dcpf_g!(ps);
            # run lsopf
            crisp_rlopf_g_orig!(ps,Pd_max);
            crisp_dcpf_g!(ps);
        else
            ps_islands = build_islands(subgraph,ps);# at some point check for changes in islands and don't run power flows if no change
            ## for every island that changed (eventually)

            for j in 1:M
                psi = ps_subset(ps,ps_islands[j]);
                    # run the dcpf
                    crisp_dcpf_g!(psi);
                    # run lsopf
                    crisp_rlopf_g_orig!(psi,Pd_max[ps_islands[j].shunt]);
                    # apply the results
                    ps.gen[ps_islands[j].gen,:Pg]  = psi.gen.Pg;
                    ps.shunt[ps_islands[j].shunt,:P] = psi.shunt.P;
                    crisp_dcpf_g!(psi);
            end
        end
        @assert abs(sum(ps.gen.Pg)-sum(ps.shunt.P))<=tolerance
        # set load shed for this time step
        load_shed[i+2] = sum(load_cost.*(Pd_max - ps.shunt[:P]));
    end
    times = [0.0;t0*1.0;Time.+t0*1.0];
    perc_load_served = (sum(load_cost.*Pd_max) .- load_shed)./sum(load_cost.*Pd_max);
    Restore = DataFrame(time = times, load_shed = load_shed, perc_load_served = perc_load_served, num_lines_out = lines_out, num_gens_out = gens_out);
    return Restore
end

function crisp_rlopf_g_orig!(ps,Pd_max)
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
        gst = (ps.gen.status .== 1);
        ng = size(ps.gen[gst,:Pg],1)
        G = bi[ps.gen[gst,:bus]]
        G_bus = sparse(G,collect(1:ng),1.,n,ng);
        Pg = ps.gen[gst,:Pg] ./ ps.baseMVA .* ps.gen[gst,:status]
        Pg_max = ps.gen[gst,:Pmax] ./ ps.baseMVA .* ps.gen[gst,:status]
        Pg_min = ps.gen[gst,:Pmin] ./ ps.baseMVA .* ps.gen[gst,:status]
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
        @variable(m1,ndPg[1:ng])
        @variable(m1,pdPg[1:ng])
        @variable(m1,ug[1:ng],Bin)
        @variable(m1,dTheta[1:n])
        # variable bounds
        @constraint(m1,-Pd.<=dPd.<=(Pd_max./ps.baseMVA - Pd));
        @constraint(m1, ug.*(Pg_min) .<= Pg+ndPg)
        @constraint(m1, ug.*(Pg_max) .>= Pg+pdPg)
        @constraint(m1, pdPg .>= 0)
        @constraint(m1, ndPg .<= 0)
        @constraint(m1,dTheta[1] == 0)
        # objective
        @objective(m1,Max,sum(dPd)+0.01*(0.1*sum(ndPg)-0.1*sum(pdPg)+sum(ug))) # serve as much load as possible
        # mapping matrix to map loads/gens to buses
        # Power balance equality constraint
        @constraint(m1,B*dTheta .== G_bus*ndPg+G_bus*pdPg-D_bus*dPd);#M_G*dPg - M_D*dPd)
        # Power flow constraints
        @constraint(m1,-flow_max .<= flow0 + Xinv.*(dTheta[F] - dTheta[T]) .<= flow_max)
        ### solve the model ###
        optimize!(m1);
        # collect/return the outputs
        sol_dPd=value.(dPd)
        sol_ndPg=value.(ndPg)
        sol_pdPg=value.(pdPg)
        dPd_star = sol_dPd.*ps.baseMVA
        dPg_star = (sol_pdPg+sol_ndPg).*ps.baseMVA
        ps.shunt.P += dPd_star; #changes ps structure
        ps.gen.Pg[gst] += dPg_star; #changes ps structure
    else
        if ((!isempty(ps.gen.status.==1) && isempty(ps.shunt)) || (isempty(ps.gen.status.==1) && !isempty(ps.shunt))
             || (isempty(ps.gen.status.==1) && isempty(ps.shunt)))
            ps.gen.Pg  .= ps.gen.Pg.*0.0;
            ps.shunt.P .= ps.shunt.P.*0.0;
        elseif !isempty(ps.gen.status.==1) && !isempty(ps.shunt)
            gst = (ps.gen.status .== 1)
            Pd = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
            Pg = ps.gen.Pg[gst] ./ ps.baseMVA .* ps.gen.status[gst]
            Pg_cap = ps.gen.Pmax[gst] ./ ps.baseMVA .* ps.gen.status[gst]
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
            ps.gen.Pg[gst]  .+= deltaPg_star;
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

## CRISP_RLSOPF_gen.jl
include("CRISP_network_gen.jl")
include("CRISP_LSOPF_gen.jl")

function RLSOPF_g!(ps,l_failures,g_failures,l_recovery_times,g_recovery_times,Pd_max;t_step = 10, t0 = 10, load_cost=0)#time is in minutes
    #add columns to keep track of the time each generator is on or off
    if sum(names(ps.gen).==:time_on) == 0
        ps.gen.time_on = zeros(length(ps.gen.Pg));
    end
    if sum(names(ps.gen).==:time_off) == 0
        ps.gen.time_off = zeros(length(ps.gen.Pg));
    end
    if sum(names(ps.gen).==:off) == 0
        ps.gen.off = zeros(length(ps.gen.Pg));
        ps.gen.off[ps.gen.Pg .== 0] .= 1;
    end
    ps.gen.time_on[ps.gen.off .== 0] .= 5000;
    # constants
    deltaT = 5; # time step in minutes;
    tolerance = 1e-6
    if sum(load_cost)==0
        load_cost = ones(length(ps.shunt.P));
    end
    l_rec_t = l_recovery_times[l_recovery_times.!=0];
    g_rec_t = g_recovery_times[g_failures.==0];
    g_rec = zeros(length(g_failures));
    g_rec[g_failures.==0] = g_rec_t;
    Time = sort([l_rec_t; g_rec_t])
    T = 0:t_step:maximum(Time)+1000;
    total_i = length(T);
    load_shed = zeros(total_i+2);
    load_shed[1] = 0;
    lines_out = zeros(total_i+2);#zeros(length(l_times)+2);
    lines_out[2] = length(l_failures) - sum(l_failures);
    gens_out = zeros(total_i+2);#zeros(length(l_times)+2);
    gens_out[2] = length(g_failures) - sum(g_failures);
    # set load shed for the step just before restoration process
    load_shed[2] = sum(load_cost.*(Pd_max - ps.shunt[:P]));
    time_gen_off = ones(length(ps.gen.Pg));
    time_gen_off[ps.gen.Pg .!=0] .=0;
    for i = 1:length(T)
        time = T[i];
        # set failed branches to status 0
        l_failures[time.>=l_recovery_times] .= 1;
        lines_out[i+2] = length(l_failures) - sum(l_failures);
        # apply to network
        ps.branch[:,:status] = l_failures;
        # set newly available generators to status 1;
        current_gen_out = ps.gen.status;
        # update generators operational
        g_failures[time.>=g_rec] .= 1;
        ps.gen.status = g_failures;
        gens_out[i+2] = length(g_failures) - sum(g_failures);
        # update time on and time off
        off = ps.gen.off .!= 0;
        on = ps.gen.off .== 0;
        ps.gen.time_off[off] .+= t_step;
        ps.gen.time_on[on] .+= t_step;
        #check for islands
        subgraph = find_subgraphs(ps);# add Int64 here hide info here
        M = Int64(findmax(subgraph)[1]);
        if M==1
            crisp_dcpf_g!(ps);
            # run lsopf
            crisp_rlopf_g!(ps,Pd_max);
            crisp_dcpf_g!(ps);
        else
            ps_islands = build_islands(subgraph,ps);# at some point check for changes in islands and don't run power flows if no change
            ## for every island that changed (eventually)

            for j in 1:M
                psi = ps_subset(ps,ps_islands[j]);
                    # run the dcpf
                    crisp_dcpf_g!(psi);
                    # run lsopf
                    crisp_rlopf_g!(psi,Pd_max[ps_islands[j].shunt]);
                    # apply the results
                    ps.gen[ps_islands[j].gen,:Pg]  = psi.gen.Pg;
                    ps.shunt[ps_islands[j].shunt,:P] = psi.shunt.P;
                    crisp_dcpf_g!(psi);
            end
        end
        @assert abs(sum(ps.gen.Pg)-sum(ps.shunt.P))<=tolerance
        # set load shed for this time step
        load_shed[i+2] = sum(load_cost.*(Pd_max - ps.shunt[:P]));
    end
    times = [0.0;t0*1.0;T.+t0*1.0];
    perc_load_served = (sum(load_cost.*Pd_max) .- load_shed)./sum(load_cost.*Pd_max);
    Restore = DataFrame(time = times, load_shed = load_shed, perc_load_served = perc_load_served, num_lines_out = lines_out, num_gens_out = gens_out);
    return Restore
end

function crisp_rlopf_g!(ps,Pd_max)
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
        gst = (ps.gen.status .== 1);
        ng = size(ps.gen[gst,:Pg],1)
        G = bi[ps.gen[gst,:bus]]
        G_bus = sparse(G,collect(1:ng),1.,n,ng);
        Pg = ps.gen[gst,:Pg] ./ ps.baseMVA .* ps.gen[gst,:status]
        Pg_max = ps.gen[gst,:Pmax] ./ ps.baseMVA .* ps.gen[gst,:status]
        Pg_min = ps.gen[gst,:Pmin] ./ ps.baseMVA .* ps.gen[gst,:status]
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
        @variable(m1,ndPg[1:ng])
        @variable(m1,pdPg[1:ng])
        @variable(m1,ug[1:ng],Bin)
        @variable(m1,dTheta[1:n])
        # variable bounds
        @constraint(m1,-Pd.<=dPd.<=(Pd_max./ps.baseMVA - Pd));
        @constraint(m1, ug.*(Pg_min) .<= Pg+ndPg)
        @constraint(m1, ug.*(Pg_max) .>= Pg+pdPg)
        @constraint(m1, pdPg .>= 0)
        @constraint(m1, ndPg .<= 0)
        @constraint(m1,dTheta[1] == 0)
        # objective
        @objective(m1,Max,sum(dPd)+0.01*(0.1*sum(ndPg)-0.1*sum(pdPg)+sum(ug))) # serve as much load as possible
        # mapping matrix to map loads/gens to buses
        # Power balance equality constraint
        @constraint(m1,B*dTheta .== G_bus*ndPg+G_bus*pdPg-D_bus*dPd);#M_G*dPg - M_D*dPd)
        # Power flow constraints
        @constraint(m1,-flow_max .<= flow0 + Xinv.*(dTheta[F] - dTheta[T]) .<= flow_max)
        ### solve the model ###
        optimize!(m1);
        # collect/return the outputs of first optimization to find which generators to turn on
        sol_dPd=value.(dPd)
        sol_ndPg=value.(ndPg)
        sol_pdPg=value.(pdPg)
        sol_ug=value.(ug)
        gen_used = Pg[sol_ug.==1];
        if sum(gen_used.==0)!=0
            ps.gen.off[gst][sol_ug.==1] .= 0;
            tst = ps.gen.time_on.*60 .>= ps.gen.StartTimeColdHr;
            gest = (gst .& tst);
            ng = size(ps.gen[gest,:Pg],1)
            G = bi[ps.gen[gest,:bus]]
            G_bus = sparse(G,collect(1:ng),1.,n,ng);
            Pg = ps.gen[gest,:Pg] ./ ps.baseMVA .* ps.gen[gest,:status]
            Pg_max = ps.gen[gest,:Pmax] ./ ps.baseMVA .* ps.gen[gest,:status]
            Pg_min = ps.gen[gest,:Pmin] ./ ps.baseMVA .* ps.gen[gest,:status]
            if any(G.<1) || any(G.>n)
                error("Bad indices in gen matrix")
            end
            m2 = Model(with_optimizer(Gurobi.Optimizer))
            # variables
            @variable(m2,dPd[1:nd])
            @variable(m2,ndPg[1:ng])
            @variable(m2,pdPg[1:ng])
            @variable(m2,ug[1:ng],Bin)
            @variable(m2,dTheta[1:n])
            # variable bounds
            @constraint(m2,-Pd.<=dPd.<=(Pd_max./ps.baseMVA - Pd));
            @constraint(m2, ug.*(Pg_min) .<= Pg+ndPg+pdPg)
            @constraint(m2, ug.*(Pg_max) .>= Pg+pdPg+ndPg)
            @constraint(m2, pdPg .>= 0)
            @constraint(m2, ndPg .<= 0)
            @constraint(m2,dTheta[1] == 0)
            # objective
            @objective(m2,Max,sum(dPd)+0.01*(0.1*sum(ndPg)-0.1*sum(pdPg)+sum(ug))) # serve as much load as possible
            # mapping matrix to map loads/gens to buses
            # Power balance equality constraint
            @constraint(m2,B*dTheta .== G_bus*ndPg+G_bus*pdPg-D_bus*dPd);#M_G*dPg - M_D*dPd)
            # Power flow constraints
            @constraint(m2,-flow_max .<= flow0 + Xinv.*(dTheta[F] - dTheta[T]) .<= flow_max)
            ### solve the model ###
            optimize!(m2);
            # collect/return the outputs
            sol_dPd=value.(dPd)
            sol_ndPg=value.(ndPg)
            sol_pdPg=value.(pdPg)
            sol_ug=value.(ug)
            dPd_star = sol_dPd.*ps.baseMVA
            dPg_star = (sol_pdPg+sol_ndPg).*ps.baseMVA
            ps.shunt.P += dPd_star; #changes ps structure
            ps.gen.Pg[gest] += dPg_star; #changes ps structure
        else
            dPd_star = sol_dPd.*ps.baseMVA
            dPg_star = (sol_pdPg+sol_ndPg).*ps.baseMVA
            ps.shunt.P += dPd_star; #changes ps structure
            ps.gen.Pg[gst] += dPg_star; #changes ps structure
        end
    else
        gst = (ps.gen[:status].==1)
        Pd = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
        nd = length(Pd)
        Pg = ps.gen.Pg[gst] ./ ps.baseMVA .* ps.gen.status[gst]
        ng = length(Pg)
        Pg_min = ps.gen.Pmin[gst] ./ ps.baseMVA .* ps.gen.status[gst]
        Pg_max = ps.gen.Pmax[gst] ./ ps.baseMVA .* ps.gen.status[gst]
        if !isempty(Pd) && !isempty(Pg)
            m = Model(with_optimizer(Gurobi.Optimizer))
            # variables
            @variable(m,dPd[1:nd])
            @variable(m,ndPg[1:ng])
            @variable(m,pdPg[1:ng])
            @variable(m,ug[1:ng],Bin)
            # variable bounds
            @constraint(m,-Pd .<= dPd .<= 0)
            @constraint(m, ug.*(Pg_min) .<= Pg+ndPg+pdPg)
            @constraint(m, ug.*(Pg_max) .>= Pg+pdPg+ndPg)
            @constraint(m, pdPg .>= 0)
            @constraint(m, ndPg .<= 0)
            # objective
            @objective(m,Max,sum(dPd)+0.01*(0.1*sum(ndPg)-0.1*sum(pdPg)+sum(ug))) # serve as much load as possible
            # mapping matrix to map loads/gens to buses
            # Power balance equality constraint
            @constraint(m,sum(Pg+ndPg+pdPg)-sum(Pd+dPd) == 0);
            ### solve the model ###
            optimize!(m);
            sol_dPd=value.(dPd)
            sol_ndPg=value.(ndPg)
            sol_pdPg=value.(pdPg)
            sol_ug=value.(ug)
            gen_used = Pg[sol_ug.==1];
            if sum(gen_used.==0)!=0
                ps.gen.off[gst][sol_ug.==1] .= 0;
                tst = ps.gen.time_on.*60 .>= ps.gen.StartTimeColdHr;
                gest = (gst .& tst);
                ng = size(ps.gen[gest,:Pg],1)
                Pg = ps.gen[gest,:Pg] ./ ps.baseMVA .* ps.gen[gest,:status]
                Pg_max = ps.gen[gest,:Pmax] ./ ps.baseMVA .* ps.gen[gest,:status]
                Pg_min = ps.gen[gest,:Pmin] ./ ps.baseMVA .* ps.gen[gest,:status]
                m = Model(with_optimizer(Gurobi.Optimizer))
                # variables
                @variable(m,dPd[1:nd])
                @variable(m,ndPg[1:ng])
                @variable(m,pdPg[1:ng])
                @variable(m,ug[1:ng],Bin)
                # variable bounds
                @constraint(m,-Pd .<= dPd .<= 0)
                @constraint(m, ug.*(Pg_min) .<= Pg+ndPg+pdPg)
                @constraint(m, ug.*(Pg_max) .>= Pg+pdPg+ndPg)
                @constraint(m, pdPg .>= 0)
                @constraint(m, ndPg .<= 0)
                # objective
                @objective(m,Max,sum(dPd)+0.01*(0.1*sum(ndPg)-0.1*sum(pdPg)+sum(ug))) # serve as much load as possible
                # mapping matrix to map loads/gens to buses
                # Power balance equality constraint
                @constraint(m,sum(Pg+ndPg+pdPg)-sum(Pd+dPd) == 0);
                ### solve the model ###
                optimize!(m);
                sol_dPd=value.(dPd)
                sol_ndPg=value.(ndPg)
                sol_pdPg=value.(pdPg)
            end
            dPd_star = sol_dPd.*ps.baseMVA
            dPg_star = (sol_pdPg+sol_ndPg).*ps.baseMVA
            ps.shunt.P .+= dPd_star;
            ps.gen.Pg[gst]  .+= dPg_star;
        else
            ps.gen.Pg  .= ps.gen.Pg.*0.0;
            ps.shunt.P .= ps.shunt.P.*0.0;
        end
    end
    #adding criteria that should produce errors if incorrect.
    @assert abs(sum(ps.shunt.P)-sum(ps.gen.Pg))<=2*tolerance
    #@assert 0.<=ps.shunt.P
    return ps
end

function RLSOPF_g_s!(ps,dt,l_failures,g_failures,l_recovery_times,g_recovery_times,Pd_max;t0 = 10, load_cost=0)#time is in minutes
    #add columns to keep track of the time each generator is on or off
    if sum(names(ps.gen).==:time_on) == 0
        ps.gen.time_on = zeros(length(ps.gen.Pg));
    end
    if sum(names(ps.gen).==:time_off) == 0
        ps.gen.time_off = zeros(length(ps.gen.Pg));
    end
    if sum(names(ps.gen).==:off) == 0
        ps.gen.off = zeros(length(ps.gen.Pg));
        ps.gen.off[ps.gen.Pg .== 0] .= 1;
    end
    ps.gen.time_on[ps.gen.off .== 0] .= 5000;
    tolerance = 1e-6
    if sum(load_cost)==0
        load_cost = ones(length(ps.shunt.P));
    end
    l_rec_t = l_recovery_times[l_recovery_times.!=0];
    g_rec_t = g_recovery_times[g_failures.==0];
    g_rec = zeros(length(g_failures));
    g_rec[g_failures.==0] = g_rec_t;
    Time = sort([l_rec_t; g_rec_t])
    T = 0:dt:maximum(Time)+1000;
    total_i = length(T);
    load_shed = zeros(total_i+2);
    load_shed[1] = 0;
    lines_out = zeros(total_i+2);#zeros(length(l_times)+2);
    lines_out[2] = length(l_failures) - sum(l_failures);
    gens_out = zeros(total_i+2);#zeros(length(l_times)+2);
    gens_out[2] = length(g_failures) - sum(g_failures);
    # set load shed for the step just before restoration process
    load_shed[2] = sum(load_cost.*(Pd_max - ps.shunt[:P]));
    time_gen_off = ones(length(ps.gen.Pg));
    time_gen_off[ps.gen.Pg .!=0] .=0;
    for i = 1:length(T)
        time = T[i];
        # set failed branches to status 0
        l_failures[time.>=l_recovery_times] .= 1;
        lines_out[i+2] = length(l_failures) - sum(l_failures);
        # apply to network
        ps.branch[:,:status] = l_failures;
        # set newly available generators to status 1;
        current_gen_out = ps.gen.status;
        # update generators operational
        g_failures[time.>=g_rec] .= 1;
        ps.gen.status = g_failures;
        gens_out[i+2] = length(g_failures) - sum(g_failures);
        # update time on and time off
        off = ps.gen.off .!= 0;
        on = ps.gen.off .== 0;
        ps.gen.time_off[off] .+= dt;
        ps.gen.time_on[on] .+= dt;
        #check for islands
        subgraph = find_subgraphs(ps);# add Int64 here hide info here
        M = Int64(findmax(subgraph)[1]);
        if M==1
            crisp_dcpf_g_s1!(ps);
            # run lsopf
            crisp_rlopf_g_s1!(ps,Pd_max,dt);
            crisp_dcpf_g_s1!(ps);
        else
            ps_islands = build_islands(subgraph,ps);# at some point check for changes in islands and don't run power flows if no change
            ## for every island that changed (eventually)

            for j in 1:M
                psi = ps_subset(ps,ps_islands[j]);
                    # run the dcpf
                    crisp_dcpf_g_s1!(psi);
                    # run lsopf
                    crisp_rlopf_g_s1!(psi,Pd_max[ps_islands[j].shunt],dt);
                    # apply the results
                    ps.gen[ps_islands[j].gen,:Pg] = psi.gen.Pg
                    ps.shunt[ps_islands[j].shunt,:P] = psi.shunt.P
                    ps.storage[ps_islands[j].storage,:E] = psi.storage.E
                    ps.storage[ps_islands[j].storage,:Ps] = psi.storage.Ps
                    crisp_dcpf_g_s1!(psi);
            end
        end
        @assert abs(sum(ps.gen.Pg)+sum(ps.storage.Ps)-sum(ps.shunt.P))<=tolerance
        # set load shed for this time step
        load_shed[i+2] = sum(load_cost.*(Pd_max - ps.shunt[:P]));
    end
    times = [0.0;t0*1.0;T.+t0*1.0];
    perc_load_served = (sum(load_cost.*Pd_max) .- load_shed)./sum(load_cost.*Pd_max);
    Restore = DataFrame(time = times, load_shed = load_shed, perc_load_served = perc_load_served, num_lines_out = lines_out, num_gens_out = gens_out);
    return Restore
end

function crisp_rlopf_g_s!(ps,Pd_max,dt)
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
        gst = (ps.gen.status .== 1);
        ng = size(ps.gen[gst,:Pg],1)
        G = bi[ps.gen[gst,:bus]]
        G_bus = sparse(G,collect(1:ng),1.,n,ng);
        Pg = ps.gen[gst,:Pg] ./ ps.baseMVA .* ps.gen[gst,:status]
        Pg_max = ps.gen[gst,:Pmax] ./ ps.baseMVA .* ps.gen[gst,:status]
        Pg_min = ps.gen[gst,:Pmin] ./ ps.baseMVA .* ps.gen[gst,:status]
        if any(G.<1) || any(G.>n)
            error("Bad indices in gen matrix")
        end
        #Pg_bus = Array(sparse(G,ones(size(G)),Pg,n,1))
        # storage data
        ns = size(ps.storage,1)
        S = bi[ps.storage[:bus]];
        S_bus = sparse(S,collect(1:ns),1.,n,ns);
        Ps = ps.storage[:Ps] ./ ps.baseMVA .* ps.storage[:status]
        E = ps.storage[:E] ./ ps.baseMVA .* ps.storage[:status]
        E_max = ps.storage[:Emax] ./ ps.baseMVA .* ps.storage[:status]
        Ps_max = ps.storage[:Psmax] ./ ps.baseMVA .* ps.storage[:status]
        Ps_min = ps.storage[:Psmin] ./ ps.baseMVA .* ps.storage[:status]
        if any(S.<1) || any(S.>n)
            error("Bad indices in sstorage matrix")
        end
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
        @variable(m1,dPs[1:ns])
        @variable(m1,ndPg[1:ng])
        @variable(m1,pdPg[1:ng])
        @variable(m1,ug[1:ng],Bin)
        @variable(m1,dTheta[1:n])
        # variable bounds
        @constraint(m1,-Pd.<=dPd.<=(Pd_max./ps.baseMVA - Pd));
        @constraint(m1, Ps_min .<= Ps+dPs .<= Ps_max)
        @constraint(m1, 0.01 .<= E + (Ps+dPs).*dt .<= E_max)
        @constraint(m1, ug.*(Pg_min) .<= Pg+ndPg)
        @constraint(m1, ug.*(Pg_max) .>= Pg+pdPg)
        @constraint(m1, pdPg .>= 0)
        @constraint(m1, ndPg .<= 0)
        @constraint(m1,dTheta[1] == 0)
        # objective
        @objective(m1,Max,sum(dPd)+0.01*(0.1*sum(ndPg)-0.1*sum(pdPg)+sum(ug))) # serve as much load as possible
        # mapping matrix to map loads/gens to buses
        # Power balance equality constraint
        @constraint(m1,B*dTheta .== G_bus*ndPg+G_bus*pdPg+S_bus*dPs-D_bus*dPd);#M_G*dPg - M_D*dPd)
        # Power flow constraints
        @constraint(m1,-flow_max .<= flow0 + Xinv.*(dTheta[F] - dTheta[T]) .<= flow_max)
        ### solve the model ###
        optimize!(m1);
        # collect/return the outputs of first optimization to find which generators to turn on
        sol_dPd=value.(dPd)
        sol_dPs=value.(dPs)
        sol_ndPg=value.(ndPg)
        sol_pdPg=value.(pdPg)
        sol_ug=value.(ug)
        gen_used = Pg[sol_ug.==1];
        if sum(gen_used.==0)!=0
            ps.gen.off[gst][sol_ug.==1] .= 0;
            tst = ps.gen.time_on.*60 .>= ps.gen.StartTimeColdHr;
            gest = (gst .& tst);
            ng = size(ps.gen[gest,:Pg],1)
            G = bi[ps.gen[gest,:bus]]
            G_bus = sparse(G,collect(1:ng),1.,n,ng);
            Pg = ps.gen[gest,:Pg] ./ ps.baseMVA .* ps.gen[gest,:status]
            Pg_max = ps.gen[gest,:Pmax] ./ ps.baseMVA .* ps.gen[gest,:status]
            Pg_min = ps.gen[gest,:Pmin] ./ ps.baseMVA .* ps.gen[gest,:status]
            if any(G.<1) || any(G.>n)
                error("Bad indices in gen matrix")
            end
            m2 = Model(with_optimizer(Gurobi.Optimizer))
            # variables
            @variable(m2,dPd[1:nd])
            @variable(m2,dPs[1:ns])
            @variable(m2,ndPg[1:ng])
            @variable(m2,pdPg[1:ng])
            @variable(m2,ug[1:ng],Bin)
            @variable(m2,dTheta[1:n])
            # variable bounds
            @constraint(m2,-Pd.<=dPd.<=(Pd_max./ps.baseMVA - Pd));
            @constraint(m2, Ps_min .<= Ps+dPs .<= Ps_max)
            @constraint(m2, 0.01 .<= E + (Ps+dPs).*dt .<= E_max)
            @constraint(m2, ug.*(Pg_min) .<= Pg+ndPg+pdPg)
            @constraint(m2, ug.*(Pg_max) .>= Pg+pdPg+ndPg)
            @constraint(m2, pdPg .>= 0)
            @constraint(m2, ndPg .<= 0)
            @constraint(m2,dTheta[1] == 0)
            # objective
            @objective(m2,Max,sum(dPd)+0.01*(0.1*sum(ndPg)-0.1*sum(pdPg)+sum(ug))) # serve as much load as possible
            # mapping matrix to map loads/gens to buses
            # Power balance equality constraint
            @constraint(m2,B*dTheta .== G_bus*ndPg+G_bus*pdPg+S_bus*dPs-D_bus*dPd);#M_G*dPg - M_D*dPd)
            # Power flow constraints
            @constraint(m2,-flow_max .<= flow0 + Xinv.*(dTheta[F] - dTheta[T]) .<= flow_max)
            ### solve the model ###
            optimize!(m2);
            # collect/return the outputs
            sol_dPd=value.(dPd)
            sol_dPs=value.(dPs)
            sol_ndPg=value.(ndPg)
            sol_pdPg=value.(pdPg)
            sol_ug=value.(ug)
            dPd_star = sol_dPd.*ps.baseMVA
            dPs_star = sol_dPs.*ps.baseMVA
            dPg_star = (sol_pdPg+sol_ndPg).*ps.baseMVA
            ps.shunt.P += dPd_star; #changes ps structure
            ps.storage.Ps += dPs_star;
            ps.storage.E += ps.storage.Ps.*dt;
            ps.gen.Pg[gest] += dPg_star; #changes ps structure
        else
            dPd_star = sol_dPd.*ps.baseMVA
            dPs_star = sol_dPs.*ps.baseMVA
            dPg_star = (sol_pdPg+sol_ndPg).*ps.baseMVA
            ps.shunt.P += dPd_star; #changes ps structure
            ps.storage.Ps += dPs_star;
            ps.storage.E += ps.storage.Ps.*dt;
            ps.gen.Pg[gst] += dPg_star; #changes ps structure
        end
    else
        gst = (ps.gen[:status].==1)
        Pd = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
        nd = length(Pd)
        Pg = ps.gen.Pg[gst] ./ ps.baseMVA .* ps.gen.status[gst]
        ng = length(Pg)
        Pg_min = ps.gen.Pmin[gst] ./ ps.baseMVA .* ps.gen.status[gst]
        Pg_max = ps.gen.Pmax[gst] ./ ps.baseMVA .* ps.gen.status[gst]
        Ps = ps.storage[:Ps] ./ ps.baseMVA .* ps.storage[:status]
        ns = length(Ps);
        E = ps.storage[:E] ./ ps.baseMVA .* ps.storage[:status]
        E_max = ps.storage[:Emax] ./ ps.baseMVA .* ps.storage[:status]
        Ps_max = ps.storage[:Ps] ./ ps.baseMVA .* ps.storage[:status]
        Ps_min = ps.storage[:Ps] ./ ps.baseMVA .* ps.storage[:status]
        if (isempty(Pd) + isempty(Pg) + isempty(Ps)) <=1
            m = Model(with_optimizer(Gurobi.Optimizer))
            # variables
            @variable(m,dPd[1:nd])
            @variable(m,dPs[1:ns])
            @variable(m,ndPg[1:ng])
            @variable(m,pdPg[1:ng])
            @variable(m,ug[1:ng],Bin)
            # variable bounds
            @constraint(m,-Pd .<= dPd .<= 0)
            @constraint(m, Ps_min .<= Ps+dPs .<= Ps_max)
            @constraint(m, 0.01 .<= E + (Ps+dPs).*dt .<= E_max)
            @constraint(m, ug.*(Pg_min) .<= Pg+ndPg+pdPg)
            @constraint(m, ug.*(Pg_max) .>= Pg+pdPg+ndPg)
            @constraint(m, pdPg .>= 0)
            @constraint(m, ndPg .<= 0)
            # objective
            @objective(m,Max,sum(dPd)+0.01*(0.1*sum(ndPg)-0.1*sum(pdPg)+sum(ug))) # serve as much load as possible
            # mapping matrix to map loads/gens to buses
            # Power balance equality constraint
            @constraint(m,sum(Pg+ndPg+pdPg)+sum(Ps+dPs)-sum(Pd+dPd) == 0);
            ### solve the model ###
            optimize!(m);
            sol_dPd=value.(dPd)
            sol_dPs=value.(dPs)
            sol_ndPg=value.(ndPg)
            sol_pdPg=value.(pdPg)
            sol_ug=value.(ug)
            gen_used = Pg[sol_ug.==1];
            if sum(gen_used.==0)!=0
                ps.gen.off[gst][sol_ug.==1] .= 0;
                tst = ps.gen.time_on.*60 .>= ps.gen.StartTimeColdHr;
                gest = (gst .& tst);
                ng = size(ps.gen[gest,:Pg],1)
                Pg = ps.gen[gest,:Pg] ./ ps.baseMVA .* ps.gen[gest,:status]
                Pg_max = ps.gen[gest,:Pmax] ./ ps.baseMVA .* ps.gen[gest,:status]
                Pg_min = ps.gen[gest,:Pmin] ./ ps.baseMVA .* ps.gen[gest,:status]
                Ps = ps.storage[:Ps] ./ ps.baseMVA .* ps.storage[:status]
                ns = length(Ps);
                E = ps.storage[:E] ./ ps.baseMVA .* ps.storage[:status]
                E_max = ps.storage[:Emax] ./ ps.baseMVA .* ps.storage[:status]
                Ps_max = ps.storage[:Ps] ./ ps.baseMVA .* ps.storage[:status]
                Ps_min = ps.storage[:Ps] ./ ps.baseMVA .* ps.storage[:status]
                m = Model(with_optimizer(Gurobi.Optimizer))
                # variables
                @variable(m,dPd[1:nd])
                @variable(m,dPs[1:ns])
                @variable(m,ndPg[1:ng])
                @variable(m,pdPg[1:ng])
                @variable(m,ug[1:ng],Bin)
                # variable bounds
                @constraint(m,-Pd .<= dPd .<= 0)
                @constraint(m, Ps_min .<= Ps+dPs .<= Ps_max)
                @constraint(m, 0.01 .<= E + (Ps+dPs).*dt .<= E_max)
                @constraint(m, ug.*(Pg_min) .<= Pg+ndPg+pdPg)
                @constraint(m, ug.*(Pg_max) .>= Pg+pdPg+ndPg)
                @constraint(m, pdPg .>= 0)
                @constraint(m, ndPg .<= 0)
                # objective
                @objective(m,Max,sum(dPd)+0.01*(0.1*sum(ndPg)-0.1*sum(pdPg)+sum(ug))) # serve as much load as possible
                # mapping matrix to map loads/gens to buses
                # Power balance equality constraint
                @constraint(m,sum(Pg+ndPg+pdPg)+sum(Ps+dPs)-sum(Pd+dPd) == 0);
                ### solve the model ###
                optimize!(m);
                sol_dPd=value.(dPd)
                sol_dPs=value.(dPs)
                sol_ndPg=value.(ndPg)
                sol_pdPg=value.(pdPg)
            end
            dPd_star = sol_dPd.*ps.baseMVA
            dPs_star = sol_dPs.*ps.baseMVA
            dPg_star = (sol_pdPg+sol_ndPg).*ps.baseMVA
            ps.shunt.P .+= dPd_star;
            ps.storage.Ps .+= dPs_star;
            ps.storage.E .+= ps.storage.Ps.*dt;
            ps.gen.Pg[gst]  .+= dPg_star;
        else
            ps.storage.Ps = ps.storage.Ps.*0.0;
            ps.gen.Pg  .= ps.gen.Pg.*0.0;
            ps.shunt.P .= ps.shunt.P.*0.0;
        end
    end
    #adding criteria that should produce errors if incorrect.
    @assert abs(sum(ps.shunt.P)-sum(ps.storage.Ps)-sum(ps.gen.Pg))<=2*tolerance
    @assert sum(ps.storage.E .<=0)==0
    #@assert 0.<=ps.shunt.P
    return ps
end

function crisp_rl_dcopf_g_s!(ps,Pd_max,dt)
    t = 1:dt:time;
    Ti = size(t,1);
    Pd_max = ps.shunt.P;
    # constants
    tolerance = 1e-6
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    if n>1
        bi = ps.bi;
        # load data
        nd = size(ps.shunt,1)
        D = bi[ps.shunt[:bus]]
        D_bus = sparse(D,collect(1:nd),1.,n,nd);
        #Pd = ps.shunt[:P] ./ ps.baseMVA .* ps.shunt[:status]
        if any(D.<1) || any(D.>n)
            error("Bad indices in shunt matrix")
        end
        #Pd_bus = Array(sparse(D,ones(size(D)),Pd,n,1))
        # gen data
        gst = (ps.gen.status .== 1);
        ng = size(ps.gen[gst,:Pg],1)
        G = bi[ps.gen[gst,:bus]]
        G_bus = sparse(G,collect(1:ng),1.,n,ng);
        #Pg = ps.gen[gst,:Pg] ./ ps.baseMVA .* ps.gen[gst,:status]
        Pg_max = ps.gen[gst,:Pmax] ./ ps.baseMVA .* ps.gen[gst,:status]
        Pg_min = ps.gen[gst,:Pmin] ./ ps.baseMVA .* ps.gen[gst,:status]
        if any(G.<1) || any(G.>n)
            error("Bad indices in gen matrix")
        end
        #Pg_bus = Array(sparse(G,ones(size(G)),Pg,n,1))
        # storage data
        ns = size(ps.storage,1)
        S = bi[ps.storage[:bus]];
        S_bus = sparse(S,collect(1:ns),1.,n,ns);
        Ps = ps.storage[:Ps] ./ ps.baseMVA .* ps.storage[:status]
        E = ps.storage[:E] ./ ps.baseMVA .* ps.storage[:status]
        E_max = ps.storage[:Emax] ./ ps.baseMVA .* ps.storage[:status]
        Ps_max = ps.storage[:Psmax] ./ ps.baseMVA .* ps.storage[:status]
        Ps_min = ps.storage[:Psmin] ./ ps.baseMVA .* ps.storage[:status]
        if any(S.<1) || any(S.>n)
            error("Bad indices in sstorage matrix")
        end
        # branch data
        brst = (ps.branch.status.==1)
        F = bi[ps.branch.f[brst]]
        T = bi[ps.branch.t[brst]]
        flow_max = ps.branch.rateA[brst]./ps.baseMVA # this could also be rateB
        Xinv = (1 ./ ps.branch.X[brst])
        B = sparse(F,T,-Xinv,n,n) +
            sparse(T,F,-Xinv,n,n) +
            sparse(T,T,+Xinv,n,n) +
            sparse(F,F,+Xinv,n,n)
        m = Model(with_optimizer(Gurobi.Optimizer))
        # variables
        @variable(m, Pd[1:nd])
        @variable(m, Ps[1:ns])
        @variable(m, Pg[1:ng])
        @variable(m, Theta[1:n])
        # variable bounds constraints
        @constraint(m, Pd .>= 0.0)
        @constraint(m, Pd .<= Pd_max)
        @constraint(m, Theta[1] .== 0);
        # power balance
        @constraint(m,B*Theta[:,k] .== G_bus*Pg-D_bus*Pd)
        # power flow
        @constraint(m, -flow_max .<= Xinv.*(Theta[F] - Theta[T]) .<= flow_max)
        # objective
        @objective(m,Max, sum(sum(Pd)))
        ## SOLVE! ##
        optimize!(m)
    else
    end
    return ps
end

## CRISP_RLOPF_gen_stor.jl
include("CRISP_network_gen.jl")
include("CRISP_LSOPF_gen_stor.jl")


function crisp_Restore(ps,l_recovery_times,g_recovery_times,dt,t_window,t0;load_cost=0)
    # constants
    tolerance = 1e-6;
    if sum(load_cost)==0
        load_cost = ones(length(ps.shunt.P));
    end
    #add columns to keep track of the time each generator is on or off
    if sum(names(ps.gen).==:time_on) == 0
        ps.gen.time_on = zeros(length(ps.gen.Pg));
    end
    if sum(names(ps.gen).==:time_off) == 0
        ps.gen.time_off = zeros(length(ps.gen.Pg));
    end
    # set time (ti)
    ti = t0;
    #save initial values
    load_shed = sum(load_cost.*(ps.shunt.P - ps.shunt.P.*ps.shunt.status));
    perc_load_served = (sum(load_cost.*ps.shunt.P) .- load_shed)./sum(load_cost.*ps.shunt.P);
    lines_out = length(ps.branch.status) - sum(ps.branch.status);
    gens_out = length(ps.gen.status) - sum(ps.gen.status);
    Restore = DataFrame(time = ti, load_shed = load_shed, perc_load_served = perc_load_served,
     lines_out = lines_out, gens_out = gens_out)
    cv = deepcopy(Restore);
    # find longest recovery
    recTime = maximum([maximum(l_recovery_times) maximum(g_recovery_times)]);
    if recTime > 1000
        dt = 60;
    end
    # set time line
    Time = t0+dt:dt:t0+recTime+dt
    for i = 1:length(Time)
        # update time
        ti = Time[i];
        # remove failures as the recovery time is reached
        ps.branch.status[ti .>= l_recovery_times] .= 1;
        ps.gen.status[ti .>= g_recovery_times] .= 1;
        # find the number of islands in ps
        subgraph = find_subgraphs(ps);# add Int64 here hide info here
        M = Int64(findmax(subgraph)[1]);
        ps_islands = build_islands(subgraph,ps);
        for j = 1:M
            psi = ps_subset(ps,ps_islands[j]);
            turn_gen_on!(psi,dt);
            crisp_lsopf_g_s!(psi,dt);
            ps.gen[ps_islands[j].gen,:Pg] = psi.gen.Pg
            ps.storage[ps_islands[j].storage,:Ps] = psi.storage.Ps
            ps.storage[ps_islands[j].storage,:E] = psi.storage.E
            ps.shunt[ps_islands[j].shunt,:status] = psi.shunt.status
            #add_changes!(ps,psi,ps_islands[j]);
        end
        # save current values
        cv.time = ti;
        cv.load_shed = sum(load_cost.*(ps.shunt.P - ps.shunt.P.*ps.shunt.status));
        cv.perc_load_served = (sum(load_cost.*ps.shunt.P) .- cv.load_shed)./sum(load_cost.*ps.shunt.P);
        cv.lines_out = length(ps.branch.status) - sum(ps.branch.status);
        cv.gens_out = length(ps.gen.status) - sum(ps.gen.status);
        append!(Restore,cv)
        println(i)
        println(cv.load_shed)
        println(ps.shunt.status)
        @assert 3*tolerance >= abs(sum(ps.shunt.P .* ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg))
    end

    # make sure that a full recovery has occued or keep iterating:
    while cv.load_shed[1] >= tolerance[1]
        # update time
        ti = ti+dt;
        # remove failures as the recovery time is reached
        ps.branch.status[ti .>= l_recovery_times] .= 1;
        ps.gen.status[ti .>= g_recovery_times] .= 1;
        # find the number of islands in ps
        subgraph = find_subgraphs(ps);# add Int64 here hide info here
        M = Int64(findmax(subgraph)[1]);
        ps_islands = build_islands(subgraph,ps);
        for j = 1:M
            psi = ps_subset(ps,ps_islands[j]);
            turn_gen_on!(psi,dt);
            crisp_lsopf_g_s!(psi,dt);
            add_changes!(ps,psi,ps_islands[j]);
        end
        @assert abs(sum(ps.shunt.P .*ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg))<=2*tolerance
        @assert sum(ps.storage.E .< 0)==0
        # save current values
        cv.time = ti;
        cv.load_shed = sum(load_cost.*(ps.shunt.P - ps.shunt.P.*ps.shunt.status));
        cv.perc_load_served = (sum(load_cost.*ps.shunt.P) .- cv.load_shed)./sum(load_cost.*ps.shunt.P);
        cv.lines_out = length(ps.branch.status) - sum(ps.branch.status);
        cv.gens_out = length(ps.gen.status) - sum(ps.gen.status);
        append!(Restore,cv)
    end
    return Restore
end

function turn_gen_on!(ps,dt)
    # constants
    tolerance = 1e-6
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    if n>1
        bi = ps.bi;
        # load data
        nd = size(ps.shunt,1)
        D = bi[ps.shunt[:bus]]
        D_bus = sparse(D,collect(1:nd),1.,n,nd);
        Pd1 = ps.shunt[:P] ./ ps.baseMVA .* ps.shunt[:status]
        if any(D.<1) || any(D.>n)
            error("Bad indices in shunt matrix")
        end
        # gen data
        gst = (ps.gen.status .== 1);
        ng = size(ps.gen[gst,:Pg],1)
        G = bi[ps.gen[gst,:bus]]
        G_bus = sparse(G,collect(1:ng),1.,n,ng);
        Pg1 = ps.gen[gst,:Pg] ./ ps.baseMVA .* ps.gen[gst,:status]
        ug1 = ones(ng); ug1[Pg1.==0] .= 0;
        #Pg = ps.gen[gst,:Pg] ./ ps.baseMVA .* ps.gen[gst,:status]
        Pg_max = ps.gen[gst,:Pmax] ./ ps.baseMVA .* ps.gen[gst,:status]
        Pg_min = ps.gen[gst,:Pmin] ./ ps.baseMVA .* ps.gen[gst,:status]
        if any(G.<1) || any(G.>n)
            error("Bad indices in gen matrix")
        end
        # storage data
        ns = size(ps.storage,1)
        S = bi[ps.storage[:bus]];
        S_bus = sparse(S,collect(1:ns),1.,n,ns);
        Ps1 = ps.storage[:Ps] ./ ps.baseMVA .* ps.storage[:status]
        E1 = ps.storage[:E] ./ ps.baseMVA .* ps.storage[:status]
        E_max = ps.storage[:Emax] ./ ps.baseMVA .* ps.storage[:status]
        Ps_max = ps.storage[:Psmax] ./ ps.baseMVA .* ps.storage[:status]
        Ps_min = ps.storage[:Psmin] ./ ps.baseMVA .* ps.storage[:status]
        if any(S.<1) || any(S.>n)
            error("Bad indices in storage matrix")
        end
        # branch data
        brst = (ps.branch.status.==1)
        F = bi[ps.branch.f[brst]]
        T = bi[ps.branch.t[brst]]
        flow_max = ps.branch.rateA[brst]./ps.baseMVA # this could also be rateB
        Xinv = (1 ./ ps.branch.X[brst])
        B = sparse(F,T,-Xinv,n,n) +
            sparse(T,F,-Xinv,n,n) +
            sparse(T,T,+Xinv,n,n) +
            sparse(F,F,+Xinv,n,n)
        m = Model(with_optimizer(Gurobi.Optimizer))
        # variables
        @variable(m, Pd[1:nd]) # demand
        @variable(m, Pg[1:ng]) # generation
        @variable(m, ug[1:ng], Bin) # generator j at time  k on == ug[j,k]=1, off == ug[j,k]=0
        @variable(m, Ps[1:ns]) # power flow into or out of storage (negative flow = charging)
        @variable(m, E[1:ns]) # energy level in battery
        @variable(m, Theta[1:n])
        # variable bounds constraints
        @constraint(m, 0.0 .<= Pd .<= ps.shunt.P./ps.baseMVA) # load served limits
        @constraint(m, Pg .<= ug.*Pg_max) # generator power limits upper
        @constraint(m, ug.*Pg_min .<= Pg) # generator power limits lower
        @constraint(m, Ps_min .<= Ps .<= Ps_max) # storage power flow
        @constraint(m, E .== (E1 + (dt/60) .*(Ps))) # storage energy at next time step
        @constraint(m, 0 .<= E .<= E_max) # storage energy
        @constraint(m, Theta[1] .== 0); # set first bus as reference bus: V angle to 0
        # power balance
        @constraint(m, B*Theta .== G_bus*Pg+S_bus*Ps-D_bus*Pd)
        # power flow limits
        @constraint(m, -flow_max .<= Xinv.*(Theta[F] - Theta[T]) .<= flow_max)
        # objective
        @objective(m, Max, 100*sum(sum(Pd)) + sum(sum(Ps)) + sum(sum(ug)));
        ## SOLVE! ##
        optimize!(m)
        sol_Pd=value.(Pd)
        sol_Ps=value.(Ps)
        sol_Pg=value.(Pg)
        sol_ug=value.(ug)
        @assert abs(sum(sol_Pd)-sum(sol_Ps)-sum(sol_Pg)) <= 3tolerance
    else
        bi = ps.bi;
        # load data
        nd = size(ps.shunt,1)
        D = bi[ps.shunt[:bus]]
        D_bus = sparse(D,collect(1:nd),1.,n,nd);
        Pd1 = ps.shunt[:P] ./ ps.baseMVA .* ps.shunt[:status]
        if any(D.<1) || any(D.>n)
            error("Bad indices in shunt matrix")
        end
        # gen data
        gst = (ps.gen.status .== 1);
        ng = size(ps.gen[gst,:Pg],1)
        G = bi[ps.gen[gst,:bus]]
        G_bus = sparse(G,collect(1:ng),1.,n,ng);
        Pg1 = ps.gen[gst,:Pg] ./ ps.baseMVA .* ps.gen[gst,:status]
        ug1 = ones(ng); ug1[Pg1.==0] .= 0;
        #Pg = ps.gen[gst,:Pg] ./ ps.baseMVA .* ps.gen[gst,:status]
        Pg_max = ps.gen[gst,:Pmax] ./ ps.baseMVA .* ps.gen[gst,:status]
        Pg_min = ps.gen[gst,:Pmin] ./ ps.baseMVA .* ps.gen[gst,:status]
        if any(G.<1) || any(G.>n)
            error("Bad indices in gen matrix")
        end
        # storage data
        ns = size(ps.storage,1)
        S = bi[ps.storage[:bus]];
        S_bus = sparse(S,collect(1:ns),1.,n,ns);
        Ps1 = ps.storage[:Ps] ./ ps.baseMVA .* ps.storage[:status]
        E1 = ps.storage[:E] ./ ps.baseMVA .* ps.storage[:status]
        E_max = ps.storage[:Emax] ./ ps.baseMVA .* ps.storage[:status]
        Ps_max = ps.storage[:Psmax] ./ ps.baseMVA .* ps.storage[:status]
        Ps_min = ps.storage[:Psmin] ./ ps.baseMVA .* ps.storage[:status]
        if any(S.<1) || any(S.>n)
            error("Bad indices in storage matrix")
        end
        m = Model(with_optimizer(Gurobi.Optimizer))
        # variables
        @variable(m, Pd[1:nd]) # demand
        @variable(m, Pg[1:ng]) # generation
        @variable(m, ug[1:ng], Bin) # generator j at time  k on == ug[j,k]=1, off == ug[j,k]=0
        @variable(m, Ps[1:ns]) # power flow into or out of storage (negative flow = charging)
        @variable(m, E[1:ns]) # energy level in battery
        # variable bounds constraints
        @constraint(m, 0.0 .<= Pd .<= ps.shunt.P./ps.baseMVA) # load served limits
        @constraint(m, Pg .<= ug.*Pg_max) # generator power limits upper
        @constraint(m, ug.*Pg_min .<= Pg) # generator power limits lower
        @constraint(m, Ps_min .<= Ps .<= Ps_max) # storage power flow
        @constraint(m, E .== (E1 + (dt/60) .*(Ps))) # storage energy at next time step
        @constraint(m, 0 .<= E .<= E_max) # storage energy
        # power balance
        @constraint(m, 0 .== G_bus*Pg+S_bus*Ps-D_bus*Pd)
        # objective
        @objective(m, Max, 100*sum(sum(Pd)) + sum(sum(ug)));
        ## SOLVE! ##
        optimize!(m)
        sol_Pd=value.(Pd)
        sol_Ps=value.(Ps)
        sol_Pg=value.(Pg)
        sol_E=value.(E)
        sol_ug=value.(ug)
        @assert abs(sum(sol_Pd)-sum(sol_Ps)-sum(sol_Pg)) <= 3tolerance
    end
    on1 =  ps.gen.time_off[gst] .>= ps.gen.minDownTimeHr[gst]
    ps.gen.time_on[gst][sol_ug.==1 .& on1] .+= dt/60;
    on2 = ps.gen.time_on[gst] .>= ps.gen.minUpTimeHr[gst]
        ps.gen.time_off[gst][on2] .= 0;
        ps.gen.Pg[gst][on2] .= ps.gen.Pmin[gst][on2]
        ps.gen.time_on[gst][on2] .= 0;
    return ps
end


## CRISP_RLOPF_movh.jl
include("CRISP_network_gen.jl")

function crisp_Restore_mh(ps,l_recovery_times,g_recovery_times,dt,t_window,t0;load_cost=0)
    # constants
    tolerance = 10^(-6);
    if sum(load_cost)==0
        load_cost = ones(length(ps.shunt.P));
    end
    ti = t0;
    #save initial values
    load_shed = sum(load_cost.*(ps.shunt.P - ps.shunt.P.*ps.shunt.status));
    perc_load_served = (sum(load_cost.*ps.shunt.P) .- load_shed)./sum(load_cost.*ps.shunt.P);
    lines_out = length(ps.branch.status) - sum(ps.branch.status);
    gens_out = length(ps.gen.status) - sum(ps.gen.status);
    Restore = DataFrame(time = ti, load_shed = load_shed, perc_load_served = perc_load_served,
     lines_out = lines_out, gens_out = gens_out)
    cv = deepcopy(Restore);
    # find longest recovery
    recTime = maximum([maximum(l_recovery_times) maximum(g_recovery_times)]);
    if recTime > 1000
        dt = 60;
        t_window = 12*60;
    end
    # set time line
    Time = t0+dt:dt:t0+recTime+dt
    for i in 1:length(Time)
        # update time
        ti = Time[i];
        # remove failures as the recovery time is reached
        ps.branch.status[ti .>= l_recovery_times] .= 1;
        ps.gen.status[ti .>= g_recovery_times] .= 1;
        # find the number of islands in ps
        subgraph = find_subgraphs(ps);# add Int64 here hide info here
        M = Int64(findmax(subgraph)[1]);
        ps_islands = build_islands(subgraph,ps);
        for j in 1:M
            psi = ps_subset(ps,ps_islands[j]);
            crisp_mh_rlopf!(psi,dt,t_window);
            ps.gen.Pg[ps_islands[j].gen] = psi.gen.Pg
            ps.gen.time_on[ps_islands[j].gen] = psi.gen.time_on
            ps.gen.time_off[ps_islands[j].gen] = psi.gen.time_off
            ps.storage.Ps[ps_islands[j].storage] = psi.storage.Ps
            ps.storage.E[ps_islands[j].storage] = psi.storage.E
            ps.shunt.status[ps_islands[j].shunt] = psi.shunt.status
            #add_changes!(ps,psi,ps_islands[j]);
        end
        # save current values
        cv.time .= ti;
        cv.load_shed .= sum(load_cost.*(ps.shunt.P - ps.shunt.P.*ps.shunt.status));
        cv.perc_load_served .= (sum(load_cost.*ps.shunt.P) .- cv.load_shed)./sum(load_cost.*ps.shunt.P);
        cv.lines_out .= length(ps.branch.status) - sum(ps.branch.status);
        cv.gens_out .= length(ps.gen.status) - sum(ps.gen.status);
        append!(Restore,cv)
        println(i)
        println(cv.load_shed)
        println(ps.shunt.status)
        println(ps.shunt.P)
        println(ps.storage.E)
        println("Pg = ")
        println(sum(ps.gen.Pg))
        println("P = ")
        println(sum(ps.shunt.P .* ps.shunt.status))
        println("Ps = ")
        println(sum(ps.storage.Ps))
        @assert 10^(-4)>=abs(sum(ps.shunt.P .* ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg))
    end

    # make sure that a full recovery has occued or keep iterating:
    while cv.load_shed[1] >= tolerance[1]
        # update time
        ti = ti+dt;
        # remove failures as the recovery time is reached
        ps.branch.status[ti .>= l_recovery_times] .= 1;
        ps.gen.status[ti .>= g_recovery_times] .= 1;
        # find the number of islands in ps
        subgraph = find_subgraphs(ps);# add Int64 here hide info here
        M = Int64(findmax(subgraph)[1]);
        ps_islands = build_islands(subgraph,ps);
        out =
        for j in 1:M
            psi = ps_subset(ps,ps_islands[j]);
            crisp_mh_rlopf!(psi,dt,t_window);
            ps.gen.Pg[ps_islands[j].gen] = psi.gen.Pg
            ps.gen.time_on[ps_islands[j].gen] = psi.gen.time_on
            ps.gen.time_off[ps_islands[j].gen] = psi.gen.time_off
            ps.storage.Ps[ps_islands[j].storage] = psi.storage.Ps
            ps.storage.E[ps_islands[j].storage] = psi.storage.E
            ps.shunt.status[ps_islands[j].shunt] = psi.shunt.status
            #add_changes!(ps,psi,ps_islands[j]);
        end
        @assert abs(sum(ps.shunt.P .*ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg))<=2*tolerance
        @assert sum(ps.storage.E .< -tolerance)==0
        # save current values
        cv.time .= ti;
        cv.load_shed .= sum(load_cost.*(ps.shunt.P - ps.shunt.P.*ps.shunt.status));
        cv.perc_load_served .= (sum(load_cost.*ps.shunt.P) .- cv.load_shed)./sum(load_cost.*ps.shunt.P);
        cv.lines_out .= length(ps.branch.status) - sum(ps.branch.status);
        cv.gens_out .= length(ps.gen.status) - sum(ps.gen.status);
        append!(Restore,cv)
    end
    return Restore
end

function crisp_mh_rlopf!(ps,dt,t_win)
    timeline = 0:dt:t_win;
    Ti = size(timeline,1);
    # constants
    tolerance = 1e-6
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    bi = ps.bi;
    # load data
    nd = size(ps.shunt,1);
    D = bi[ps.shunt.bus];
    D_bus = sparse(D,collect(1:nd),1.,n,nd);
    Pdmax = (ps.shunt.P ./ ps.baseMVA).* ps.shunt.status;
    Pd1 = (ps.shunt.P ./ ps.baseMVA) .* ps.shunt.status;
    if any(D.<1) || any(D.>n)
        error("Bad indices in shunt matrix")
    end
    # gen data
    gst = (ps.gen.status .== 1);
    ng = size(ps.gen.Pg[gst],1)
    G = bi[ps.gen.bus[gst]]
    G_bus = sparse(G,collect(1:ng),1.,n,ng);
    Pg1 = (ps.gen.Pg[gst] ./ ps.baseMVA) .* ps.gen.status[gst]
    ug1 = falses(ng); ug1[Pg1.!=0] .= true;
    Pg_max = (ps.gen.Pmax[gst] ./ ps.baseMVA) .* ps.gen.status[gst]
    Pg_min = (ps.gen.Pmin[gst] ./ ps.baseMVA) .* ps.gen.status[gst]
    RR = ps.gen.RampRateMWMin[gst] .* dt ./ ps.baseMVA .* ps.gen.status[gst]
    T_SU =  Int64.(round.(ps.gen.minUpTimeHr[gst] .*60 ./dt .* ps.gen.status[gst])) # number of time steps to turn on
    T_SD =  Int64.(round.(ps.gen.minDownTimeHr[gst] .*60 ./dt ./ ps.baseMVA .* ps.gen.status[gst])) # number of time steps to turn off
    t_off = Int64.(round.(ps.gen.time_off[gst]./dt ./ ps.baseMVA .* ps.gen.status[gst]))
    t_on = Int64.(round.(ps.gen.time_on[gst]./dt ./ ps.baseMVA .* ps.gen.status[gst]))
    if any(G.<1) || any(G.>n)
        error("Bad indices in gen matrix")
    end
    # storage data
    ns = size(ps.storage,1)
    S = bi[ps.storage.bus];
    S_bus = sparse(S,collect(1:ns),1.,n,ns);
    Ps1 = (ps.storage.Ps ./ ps.baseMVA) .* ps.storage.status
    E1 = (ps.storage.E ./ ps.baseMVA) .* ps.storage.status
    E_max = (ps.storage.Emax ./ ps.baseMVA) .* ps.storage.status
    Ps_max = (ps.storage.Psmax ./ ps.baseMVA) .* ps.storage.status
    Ps_min = (ps.storage.Psmin ./ ps.baseMVA) .* ps.storage.status
    if any(S.<1) || any(S.>n)
        error("Bad indices in storage matrix")
    end
    if n>1
        # branch data
        brst = (ps.branch.status.==1)
        F = bi[ps.branch.f[brst]]
        T = bi[ps.branch.t[brst]]
        flow_max = ps.branch.rateA[brst]./ps.baseMVA # this could also be rateB
        Xinv = (1 ./ ps.branch.X[brst])
        B = sparse(F,T,-Xinv,n,n) +
            sparse(T,F,-Xinv,n,n) +
            sparse(T,T,+Xinv,n,n) +
            sparse(F,F,+Xinv,n,n);
        # vector that depreciates the value of later elements in objective
        C_time = exp.(0:-1:(-Ti+1))';
        m = Model(with_optimizer(Gurobi.Optimizer))
        # variables
        @variable(m, Pd[1:nd, 1:Ti]) # demand
        @variable(m, Pg[1:ng, 1:Ti]) # generation
        @variable(m, ug[1:ng, 1:Ti], Bin) # generator j at time  k on == ug[j,k]=1, off == ug[j,k]=0
        @variable(m, gon[1:ng, 1:Ti], Bin) # generator j at time  k on == ug[j,k]=1, off == ug[j,k]=0
        @variable(m, goff[1:ng, 1:Ti], Bin) # generator j at time  k on == ug[j,k]=1, off == ug[j,k]=0
        @variable(m, Ps[1:ns, 1:Ti]) # power flow into or out of storage (negative flow = charging)
        @variable(m, E[1:ns, 1:Ti]) # energy level in battery
        @variable(m, Theta[1:n, 1:Ti])
        # first time step constraints
        for k=1
            for d in 1:nd
                fix(Pd[d,k], Pd1[d], force = true) #Pd1[d]
            end
            for s in 1:ns
                fix(Ps[s,k], Ps1[s], force = true)
                fix(E[s,k], E1[s], force = true)
            end
            for g in 1:ng
                fix(Pg[g,k], Pg1[g], force = true) #Pg1[g]
                @constraint(m, ug[g,k] .== ug1[g]) #ug1[g]
            end
        end
        # variable bounds constraints
        @constraint(m, stPdcon[k=2:Ti], 0.0 .<= Pd[:,k] .<= Pdmax) # load served limits
        @constraint(m, stPscon[k=2:Ti], Ps_min .<= Ps[:,k] .<= Ps_max) # storage power flow
        @constraint(m, stEPscon[k=2:Ti], E[:,k] .== (E[:,k-1] + ((dt/60) .*(Ps[:,k])))) # storage energy at next time step
        @constraint(m, stEcon[k=2:Ti], 0 .<= (E[:,k]) .<= E_max) # storage energy
        @constraint(m, genPgucon[k=2:Ti], Pg[:,k] .<= ug[:,k].*Pg_max) # generator power limits upper
        @constraint(m, genPglcon[k=2:Ti], ug[:,k].*Pg_min .<= Pg[:,k]) # generator power limits lower
        @constraint(m, genOnOffcon[k=2:Ti], ug[:,k] .<= ug[:,k-1] + gon[:,k-1] - goff[:,k-1]) # generator on and off constraint
        #power balance
        @constraint(m, PBcon[k=2:Ti], B*Theta[:,k] .== G_bus*Pg[:,k]+S_bus*Ps[:,k]-D_bus*Pd[:,k])
        # power flow limits
        @constraint(m, PFcon[k=2:Ti], -flow_max .<= Xinv .* (Theta[F,k] - Theta[T,k]) .<= flow_max)
        #Sum_Pd = 0;
        #Sum_ug = 0;
        for g in 1:ng
            for k in 2:Ti
                # Start-up time constraint
                #Summation term in the start-up time constraint:
                sum_ug_on = 0;
                for i in (k-(T_SU[g]+T_SD[g])):k # should this be (T_SU[g]+T_SD[g])
                    if i <= 0
                        sum_ug = t_off[g]+t_on[g];
                    elseif ug[g,i].==0
                        sum_ug +=1;
                    else

                        #fix(gon[g,k], 0,force = true);
                    end
                end
                if sum_ug_on >= T_SU[g]+T_SD[g]
                    @constraint(m, gon[g,k] == 0);
                end
                #Shut-down time constraint
                sum_ug_off = 0;
                for l in k:k+T_SD[g]
                    #Summation term in the shutdown time constraint
                    if l > Ti
                        sum_ug_off += 1;
                    elseif ug[g,l] == 0
                        sum_ug_off += 1;
                    else
                        @constraint(m, goff[g,k] == 0);
                        #fix(goff[g,k], 0,force = t);
                    end
                end
                #ramp rate constraints
                if ((goff[g,k-1] == 0) & (gon[g,k-1] == 0))
                    @constraint(m, -RR .<= Pg[g,k-1]-Pg[g,k] .<= RR)
                end
            end
        end
        #@constraint(m, [g=1:ng, k=1:Ti],   sum(1 .- ug[g,k-T_SU[g]:k]) .>= T_SU[g].*gon[g,k]) # generator power start up
        #@constraint(m, [g=1:ng, k=1:Ti],   sum(1 .- ug[g,k:k+T_SD[g]]) .>= T_SD[g].*goff[g,k]) # generator power shut down
        @constraint(m, Theta[1,:] .== 0); # set first bus as reference bus: V angle to 0
        # objective
        @objective(m, Max, 100*sum(Pd*C_time') + sum(ug*C_time'));
        ## SOLVE! ##
        optimize!(m)
        sol_Pd=value.(Pd)[:,2]
        sol_Ps=value.(Ps)[:,2]
        sol_Pg=value.(Pg)[:,2]
        sol_E=value.(E)[:,2]
        sol_ug=value.(ug)
        sol_gon = value.(gon)
        sol_goff = value.(goff)
        @assert abs(sum(ps.shunt.P.*ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg[gst]))<=2*tolerance
        @assert sum(ps.storage.E .< -tolerance)==0
        dPd_star = (Vector(sol_Pd).*ps.baseMVA)./ps.shunt.P # % load served
        dPs_star = Vector(sol_Ps).*ps.baseMVA
        dPg_star = Vector(sol_Pg).*ps.baseMVA
        dE_star = Vector(sol_E).*ps.baseMVA
    else
        # vector that depreciates the value of later elements in objective
        C_time = exp.(0:-1:(-Ti+1))';
        m = Model(with_optimizer(Gurobi.Optimizer))
        # variables
        @variable(m, Pd[1:nd, 1:Ti]) # demand
        @variable(m, Pg[1:ng, 1:Ti]) # generation
        @variable(m, ug[1:ng, 1:Ti], Bin) # generator j at time  k on == ug[j,k]=1, off == ug[j,k]=0
        @variable(m, gon[1:ng, 1:Ti], Bin) # generator j at time  k on == ug[j,k]=1, off == ug[j,k]=0
        @variable(m, goff[1:ng, 1:Ti], Bin) # generator j at time  k on == ug[j,k]=1, off == ug[j,k]=0
        @variable(m, Ps[1:ns, 1:Ti]) # power flow into or out of storage (negative flow = charging)
        @variable(m, E[1:ns, 1:Ti]) # energy level in battery
        # first time step constraints
        for k=1
            for d in 1:nd
                fix(Pd[d,k], 0.0, force = true) #Pd1[d]
            end
            for s in 1:ns
                fix(Ps[s,k], Ps1[s], force = true)
                fix(E[s,k], E1[s], force = true)
            end
            for g in 1:ng
                fix(Pg[g,k], 0.0, force = true) #Pd1[d]
                #fix(ug[g,k], ug1[g], force = true)
                @constraint(m, ug[g,k] .== false) #ug1[g]
                #@constraint(m, Pg[g,1] .== Pg1)
            end
        end
        # variable bounds constraints
        @constraint(m, stPdcon[k=2:Ti], 0.0 .<= Pd[:,k] .<= Pdmax) # load served limits
        @constraint(m, stPscon[k=2:Ti], Ps_min .<= Ps[:,k] .<= Ps_max) # storage power flow
        @constraint(m, stEPscon[k=2:Ti], E[:,k] .== (E[:,k-1] + ((dt/60) .*(Ps[:,k])))) # storage energy at next time step
        @constraint(m, stEcon[k=2:Ti], 0.01 .<= (E[:,k]) .<= E_max) # storage energy
        @constraint(m, genPgucon[k=2:Ti], Pg[:,k] .<= ug[:,k].*Pg_max) # generator power limits upper
        @constraint(m, genPglcon[k=2:Ti], ug[:,k].*Pg_min .<= Pg[:,k]) # generator power limits lower
        @constraint(m, genOnOffcon[k=2:Ti], ug[:,k] .<= ug[:,k-1] + gon[:,k-1] - goff[:,k-1]) # generator on and off constraint
        #power balance
        @constraint(m, PBcon[k=2:Ti], 0.0 .== G_bus*Pg[:,k]+S_bus*Ps[:,k]-D_bus*Pd[:,k])
        for g in 1:ng
            for k in 2:Ti
                # Start-up time constraint
                #Summation term in the start-up time constraint:
                sum_ug_on = 0;
                for i in k-(T_SU[g]):k # should this be (T_SU[g]+T_SD[g])
                    if i <= 0
                        sum_ug = t_on[g];
                    elseif ug[g,i].==0
                        sum_ug +=1;
                    else
                        @constraint(m, gon[g,k] == 0);
                    end
                end
                #Shut-down time constraint
                sum_ug_off = 0;
                for l in k:k+T_SD[g]
                    #Summation term in the shutdown time constraint
                    if l > Ti
                        sum_ug_off += 1;
                    elseif ug[g,l] == 0
                        sum_ug_off += 1;
                    else
                        @constraint(m, goff[g,k] == 0);
                    end
                end
                #ramp rate constraints
                if ((goff[g,k-1] == 0) & (gon[g,k-1] == 0))
                    @constraint(m, -RR .<= Pg[g,k-1]-Pg[g,k] .<= RR)
                end
            end
        end
        @constraint(m, Theta[1,:] .== 0); # set first bus as reference bus: V angle to 0
        # objective
        println(size(Pd))
        println(size(C_time))
        println(size(ug))
        @objective(m, Max, 100*sum(Pd*C_time') + sum(ug*C_time'));
        ## SOLVE! ##
        optimize!(m)
        sol_Pd=value.(Pd)[:,2]
        sol_Ps=value.(Ps)[:,2]
        sol_Pg=value.(Pg)[:,2]
        sol_E=value.(E)[:,2]
        sol_ug=value.(ug)
        sol_gon = value.(gon)
        sol_goff = value.(goff)
        @assert abs(sum(ps.shunt.P.*ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg[gst]))<=2*tolerance
        @assert sum(ps.storage.E .< -tolerance)==0
        dPd_star = (Vector(sol_Pd).*ps.baseMVA)./ps.shunt.P # % load served
        dPs_star = Vector(sol_Ps).*ps.baseMVA
        dPg_star = Vector(sol_Pg).*ps.baseMVA
        dE_star = Vector(sol_E).*ps.baseMVA
    end
    # add changes ps/psi structure
    ps.shunt.status = dPd_star;
    ps.storage.Ps = dPs_star;
    ps.storage.E = dE_star;
    ps.gen.Pg[gst] = dPg_star;
    # find turn on and off time
    ps.gen.time_off[ps.gen.Pg .== 0] .+= dt;
    ps.gen.time_off[ps.gen.Pg .!= 0] .= 0.0;
    return ps
end
