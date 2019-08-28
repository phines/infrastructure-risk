# moving horizon model of restoration process after blackout
# includes generator warm up and shut down times and storage
# still needs variable weather and load data to have appropriate hourly responses
using DataFrames;
using JuMP
#using Gurobi
using Cbc;
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
            turn_gen_on!(psi,dt);
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
    g1 = (ps.gen.status .== 1);
    g2 = (ps.gen.Pg .!=0);
    g3 = (ps.gen.time_on .!= 0)
    gst = g1 .| g2;
    gf = .!gst .& g3
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
        #m = Model(with_optimizer(Gurobi.Optimizer))
	m = Model(with_optimizer(Cbc.Optimizer))
        # variables
        @variable(m, Pd[1:nd, 1:Ti]) # demand
        @variable(m, Pg[1:ng, 1:Ti]) # generation
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
                fix(Pg[g,k], Pg1[g], force = true) #Pg1[gst]
            end
        end
        # variable bounds constraints
        @constraint(m, stPdcon[k=2:Ti], 0.0 .<= Pd[:,k] .<= Pdmax) # load served limits
        @constraint(m, stPscon[k=2:Ti], Ps_min .<= Ps[:,k] .<= Ps_max) # storage power flow
        @constraint(m, stEPscon[k=2:Ti], E[:,k] .== (E[:,k-1] + ((dt/60) .*(Ps[:,k])))) # storage energy at next time step
        @constraint(m, stEcon[k=2:Ti], 0 .<= (E[:,k]) .<= E_max) # storage energy
        @constraint(m, genPgucon[k=2:Ti], Pg[:,k] .<= Pg_max) # generator power limits upper
        @constraint(m, genPglcon[k=2:Ti], 0 .<= Pg[:,k]) # generator power limits lower
        #power balance
        @constraint(m, PBcon[k=2:Ti], B*Theta[:,k] .== G_bus*Pg[:,k]+S_bus*Ps[:,k]-D_bus*Pd[:,k])
        # power flow limits
        @constraint(m, PFcon[k=2:Ti], -flow_max .<= Xinv .* (Theta[F,k] - Theta[T,k]) .<= flow_max)
        @constraint(m, Theta[1,:] .== 0); # set first bus as reference bus: V angle to 0
        # objective
        @objective(m, Max, 100*sum(Pd*C_time'));
        ## SOLVE! ##
        optimize!(m)
        sol_Pd=value.(Pd)[:,2]
        sol_Ps=value.(Ps)[:,2]
        sol_Pg=value.(Pg)[:,2]
        sol_E=value.(E)[:,2]
        @assert abs(sum(ps.shunt.P.*ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg[gst]))<=2*tolerance
        @assert sum(ps.storage.E .< -tolerance)==0
        dPd_star = (Vector(sol_Pd).*ps.baseMVA)./ps.shunt.P # % load served
        dPs_star = Vector(sol_Ps).*ps.baseMVA
        dPg_star = Vector(sol_Pg).*ps.baseMVA
        dE_star = Vector(sol_E).*ps.baseMVA
    else
        # vector that depreciates the value of later elements in objective
        C_time = exp.(0:-1:(-Ti+1))';
        #m = Model(with_optimizer(Gurobi.Optimizer))
	m = Model(with_optimizer(Cbc.Optimizer))
        # variables
        @variable(m, Pd[1:nd, 1:Ti]) # demand
        @variable(m, Pg[1:ng, 1:Ti]) # generation
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
                fix(Pg[gst,k], 0.0, force = true) #Pd1[d]
                #fix(ug[gst,k], ug1[gst], force = true)
                @constraint(m, ug[gst,k] .== ug1[gst])
                #@constraint(m, Pg[gst,1] .== Pg1)
            end
        end
        # variable bounds constraints
        @constraint(m, stPdcon[k=2:Ti], 0.0 .<= Pd[:,k] .<= Pdmax) # load served limits
        @constraint(m, stPscon[k=2:Ti], Ps_min .<= Ps[:,k] .<= Ps_max) # storage power flow
        @constraint(m, stEPscon[k=2:Ti], E[:,k] .== (E[:,k-1] + ((dt/60) .*(Ps[:,k])))) # storage energy at next time step
        @constraint(m, stEcon[k=2:Ti], 0.01 .<= (E[:,k]) .<= E_max) # storage energy
        @constraint(m, genPgucon[k=2:Ti], Pg[:,k] .<= ug[:,k].*Pg_max) # generator power limits upper
        @constraint(m, genPglcon[k=2:Ti], 0 .<= Pg[:,k]) # generator power limits lower
        #power balance
        @constraint(m, PBcon[k=2:Ti], 0.0 .== G_bus*Pg[:,k]+S_bus*Ps[:,k]-D_bus*Pd[:,k])
        #
        @constraint(m, Theta[1,:] .== 0); # set first bus as reference bus: V angle to 0
        # objective
        println(size(Pd))
        println(size(C_time))
        println(size(ug))
        @objective(m, Max, 100*sum(Pd*C_time'));
        ## SOLVE! ##
        optimize!(m)
        sol_Pd=value.(Pd)[:,2]
        sol_Ps=value.(Ps)[:,2]
        sol_Pg=value.(Pg)[:,2]
        sol_E=value.(E)[:,2]
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
        gs = (ps.gen.Pg .!= 0);
        gst = (ps.gen.status .== 1);
        g = gs .& gst;
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
        gs = (ps.gen.Ps .!= 0);
        gst = (ps.gen.status .== 1);
        g = gs .& gst;
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
