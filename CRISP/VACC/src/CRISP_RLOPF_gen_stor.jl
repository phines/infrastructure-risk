# moving horizon model of restoration process after blackout
# includes generator warm up and shut down times and storage
# still needs variable weather and load data to have appropriate hourly responses
using DataFrames;
using JuMP
#using Gurobi
using Cbc;
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
        #m = Model(with_optimizer(Gurobi.Optimizer))
	m = Model(with_optimizer(Cbc.Optimizer))
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
        @constraint(m, E .== (E1 - (dt/60) .*(Ps))) # storage energy at next time step
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
        #m = Model(with_optimizer(Gurobi.Optimizer))
	m = Model(with_optimizer(Cbc.Optimizer))
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
        @constraint(m, E .== (E1 - (dt/60) .*(Ps))) # storage energy at next time step
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
