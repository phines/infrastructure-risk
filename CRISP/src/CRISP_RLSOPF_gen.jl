using DataFrames;
using SpecialFunctions;
using JuMP;
using Clp;
using Gurobi;
using Cbc;
include("CRISP_LSOPF_gen.jl")
include("CRISP_network_gen.jl")

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
            crisp_dcpf_g_s!(ps);
            # run lsopf
            crisp_rlopf_g_s!(ps,Pd_max,dt);
            crisp_dcpf_g_s!(ps);
        else
            ps_islands = build_islands(subgraph,ps);# at some point check for changes in islands and don't run power flows if no change
            ## for every island that changed (eventually)

            for j in 1:M
                psi = ps_subset(ps,ps_islands[j]);
                    # run the dcpf
                    crisp_dcpf_g_s!(psi);
                    # run lsopf
                    crisp_rlopf_g_s!(psi,Pd_max[ps_islands[j].shunt],dt);
                    # apply the results
                    ps.gen[ps_islands[j].gen,:Pg] = psi.gen.Pg
                    ps.shunt[ps_islands[j].shunt,:P] = psi.shunt.P
                    ps.storage[ps_islands[j].storage,:E] = psi.storage.E
                    ps.storage[ps_islands[j].storage,:Ps] = psi.storage.Ps
                    crisp_dcpf_g_s!(psi);
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
        Ps_max = ps.storage[:Ps] ./ ps.baseMVA .* ps.storage[:status]
        Ps_min = ps.storage[:Ps] ./ ps.baseMVA .* ps.storage[:status]
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
            ns = size(ps.storage,1)
            S = bi[ps.storage[:bus]];
            S_bus = sparse(S,collect(1:ns),1.,n,ns);
            Ps = ps.storage[:Ps] ./ ps.baseMVA .* ps.storage[:status]
            E = ps.storage[:E] ./ ps.baseMVA .* ps.storage[:status]
            E_max = ps.storage[:Emax] ./ ps.baseMVA .* ps.storage[:status]
            Ps_max = ps.storage[:Ps] ./ ps.baseMVA .* ps.storage[:status]
            Ps_min = ps.storage[:Ps] ./ ps.baseMVA .* ps.storage[:status]
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
