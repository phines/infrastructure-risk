# moving horizon model of restoration process after blackout
# includes generator warm up and shut down times and storage
# still needs variable weather and load data to have appropriate hourly responses
using DataFrames;
using JuMP
using Gurobi
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
            crisp_mh_rlopf!(psi,dt,t_window);
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

function crisp_mh_rlopf1!(ps,dt,t_win)
    t = 0:dt:t_win;
    Ti = size(t,1);
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
        Pd1 = (ps.shunt[:P] ./ ps.baseMVA) .* ps.shunt[:status]
        if any(D.<1) || any(D.>n)
            error("Bad indices in shunt matrix")
        end
        # gen data
        gst = (ps.gen.status .== 1);
        ng = size(ps.gen[gst,:Pg],1)
        G = bi[ps.gen[gst,:bus]]
        G_bus = sparse(G,collect(1:ng),1.,n,ng);
        Pg1 = (ps.gen[gst,:Pg] ./ ps.baseMVA) .* ps.gen[gst,:status]
        ug1 = ones(ng); ug1[Pg1.==0] .= 0;
        Pg_max = (ps.gen[gst,:Pmax] ./ ps.baseMVA) .* ps.gen[gst,:status]
        Pg_min = (ps.gen[gst,:Pmin] ./ ps.baseMVA) .* ps.gen[gst,:status]
        RR = ps.gen[gst,:RampRateMWMin] .* dt ./ ps.baseMVA .* ps.gen[gst,:status]
        T_SU =  ps.gen[gst,:minUpTimeHr] .*60 ./dt .* ps.gen[gst,:status] # number of time steps to turn on
        T_SD =  ps.gen[gst,:minDownTimeHr] .*60 ./dt ./ ps.baseMVA .* ps.gen[gst,:status] # number of time steps to turn off
        Tp = maximum(T_SU);
        if any(G.<1) || any(G.>n)
            error("Bad indices in gen matrix")
        end
        # storage data
        ns = size(ps.storage,1)
        S = bi[ps.storage[:bus]];
        S_bus = sparse(S,collect(1:ns),1.,n,ns);
        Ps1 = (ps.storage[:Ps] ./ ps.baseMVA) .* ps.storage[:status]
        E1 = (ps.storage[:E] ./ ps.baseMVA) .* ps.storage[:status]
        E_max = (ps.storage[:Emax] ./ ps.baseMVA) .* ps.storage[:status]
        Ps_max = (ps.storage[:Psmax] ./ ps.baseMVA) .* ps.storage[:status]
        Ps_min = (ps.storage[:Psmin] ./ ps.baseMVA) .* ps.storage[:status]
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
        @variable(m, Pd[1:nd, 1:Ti]) # demand
        @variable(m, Pg[1:ng, 1:Ti]) # generation
        @variable(m, ug[1:ng, -Tp:Ti], Bin) # generator j at time  k on == ug[j,k]=1, off == ug[j,k]=0
        @variable(m, gon[1:ng, -Tp:Ti], Bin) # generator j at time  k on == ug[j,k]=1, off == ug[j,k]=0
        @variable(m, goff[1:ng, -Tp:Ti], Bin) # generator j at time  k on == ug[j,k]=1, off == ug[j,k]=0
        @variable(m, Ps[1:ns, 1:Ti]) # power flow into or out of storage (negative flow = charging)
        @variable(m, E[1:ns, 1:Ti]) # energy level in battery
        @variable(m, Theta[1:n, 1:Ti])
        # variable bounds constraints
        @constraint(m, [k=1:Ti], 0.0 .<= Pd[:,k] .<= (ps.shunt.P./ps.baseMVA)) # load served limits
        @constraint(m, [k=1:Ti], Pg[:,k] .<= ug[:,k].*Pg_max) # generator power limits upper
        @constraint(m, [k=1:Ti], ug[:,k].*Pg_min .<= Pg[:,k]) # generator power limits lower
        #@constraint(m, [k=1:Ti], ug[:,k] .<= ug[:,k-1] + gon[:,k-1] - goff[:,k-1]) # generator on and off constraint
        # Shutdown time constraint
        for g in 1:ng
            for k in 1:Ti
                #Summation term in the shutdown time constraint:
                if sum(1 .- ug[g,k-T_SU[g]:k]) .>= T_SU[g])
                else
                    @constraint(m, gon[g,k] == 0)
                end
                #=
                for idx in k-ps.gen.ShutDownTime[g]:k-1
                    if idx>0
                        gen_state_exp = gen_state_exp + (1-ug[idx,g])
                        counter+=1
                    elseif idx<=0 && size(prev_run_data) > (0,0)  && size(prev_run_data,1) > abs(idx)
                        gen_state_exp = gen_state_exp + (1-prev_run_data[Symbol("Ug$g")][end+idx])
                        counter+=1
                    else #idx<=0 && size(prev_run_data) == (0,0)
                        #Ignore previous generator states
                    end
                end
                #counter == ps.gen.ShutDownTime[g] under most conditions, but not all
                #NaNs don't affect Cbc, but Gurobi gets angry
                if counter != 0 || gen_state_exp != 0
                    @constraint(uc, us[k,g]<=1/counter*gen_state_exp)
                end =#
            end
        end
        #@constraint(m, [g=1:ng, k=1:Ti],   sum(1 .- ug[g,k-T_SU[g]:k]) .>= T_SU[g].*gon[g,k]) # generator power start up
        #@constraint(m, [g=1:ng, k=1:Ti],   sum(1 .- ug[g,k:k+T_SD[g]]) .>= T_SU[g].*goff[g,k]) # generator power shut down
        @constraint(m, [k=1:Ti], Ps_min .<= Ps[:,k] .<= Ps_max) # storage power flow
        @constraint(m, [j=2:Ti], E[:,j] .== (E[:,j-1] + ((dt/60) .*(Ps[:,j])))) # storage energy at next time step
        @constraint(m, [i=1:Ti], 0 .<= (E[:,i]) .<= E_max) # storage energy
        @constraint(m, Theta[1,:] .== 0); # set first bus as reference bus: V angle to 0
        # set starting point (time at step 0 == k=1);
        @constraint(m, Pd[:,1] .== Pd1)
        @constraint(m, ug[:,1] .== ug1)
        @constraint(m, gon[:,-Tp:1] .== gon1)
        @constraint(m, goff[:,-Tp:1] .== goff1)
        @constraint(m, Pg[:,1] .== Pg1)
        @constraint(m, Ps[:,1] .== Ps1)
        @constraint(m, E[:,1] .== E1)
        # power balance
        @constraint(m, [k=1:Ti], B*Theta[:,k] .== G_bus*Pg[:,k]+S_bus*Ps[:,k]-D_bus*Pd[:,k])
        # power flow limits
        #@constraint(m, [k=1:Ti], -flow_max .<= Xinv.*(Theta[F,k] - Theta[T,k]) .<= flow_max)
        C_time = exp.(-1:-1:-Ti); # depreciates value of
        # objective
        @objective(m, Max, 100*sum(Pd*C_time) + sum(ug*C_time);
        ## SOLVE! ##
        optimize!(m)
        sol_Pd=value.(Pd)[:,2]
        sol_Ps=value.(Ps)[:,2]
        sol_Pg=value.(Pg)[:,2]
        sol_E=value.(E)[:,2]
        sol_ug=value.(ug)[:,2]
        @assert abs(sum(ps.shunt.P.*ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg[gst]))<=2*tolerance
        @assert sum(ps.storage.E .< 0)==0
        dPd_star = (sol_Pd.*ps.baseMVA)./ps.shunt.P # % load served
        dPs_star = sol_Ps.*ps.baseMVA
        dPg_star = sol_Pg.*ps.baseMVA
        dE_star = sol_E.*ps.baseMVA
    else
        bi = ps.bi;
        # load data
        nd = size(ps.shunt,1)
        D = bi[ps.shunt[:bus]]
        D_bus = sparse(D,collect(1:nd),1.,n,nd);
        Pd1 = (ps.shunt[:P] ./ ps.baseMVA) .* ps.shunt[:status]
        if any(D.<1) || any(D.>n)
            error("Bad indices in shunt matrix")
        end
        # gen data
        gst = (ps.gen.status .== 1);
        ng = size(ps.gen[gst,:Pg],1)
        G = bi[ps.gen[gst,:bus]]
        G_bus = sparse(G,collect(1:ng),1.,n,ng);
        Pg1 = (ps.gen[gst,:Pg] ./ ps.baseMVA) .* ps.gen[gst,:status]
        ug1 = ones(ng); ug1[Pg1.==0] .= 0;
        RR = ps.gen[gst,:RampRateMWMin] .* dt ./ ps.baseMVA .* ps.gen[gst,:status]
        Pg_max = (ps.gen[gst,:Pmax] ./ ps.baseMVA) .* ps.gen[gst,:status]
        Pg_min = (ps.gen[gst,:Pmin] ./ ps.baseMVA) .* ps.gen[gst,:status]
        T_SU =  ps.gen[gst,:minUpTimeHr] ./ ps.baseMVA .* ps.gen[gst,:status]
        T_SD =  ps.gen[gst,:minDownTimeHr] ./ ps.baseMVA .* ps.gen[gst,:status]
        if any(G.<1) || any(G.>n)
            error("Bad indices in gen matrix")
        end
        # storage data
        ns = size(ps.storage,1)
        S = bi[ps.storage[:bus]];
        S_bus = sparse(S,collect(1:ns),1.,n,ns);
        Ps1 = (ps.storage[:Ps] ./ ps.baseMVA) .* ps.storage[:status]
        E1 = (ps.storage[:E] ./ ps.baseMVA) .* ps.storage[:status]
        E_max = (ps.storage[:Emax] ./ ps.baseMVA) .* ps.storage[:status]
        Ps_max = (ps.storage[:Psmax] ./ ps.baseMVA) .* ps.storage[:status]
        Ps_min = (ps.storage[:Psmin] ./ ps.baseMVA) .* ps.storage[:status]
        if any(S.<1) || any(S.>n)
            error("Bad indices in storage matrix")
        end
        m = Model(with_optimizer(Gurobi.Optimizer))
        # variables
        @variable(m, Pd[1:nd, 1:Ti]) # demand
        @variable(m, Pg[1:ng, 1:Ti]) # generation
        @variable(m, ug[1:ng, 1:Ti], Bin) # generator j at time  k on == ug[j,k]=1, off == ug[j,k]=0
        @variable(m, Ps[1:ns, 1:Ti]) # power flow into or out of storage (negative flow = charging)
        @variable(m, E[1:ns,1:Ti+1]) # energy level in battery
        # variable bounds constraints
        @constraint(m, [k=1:Ti], 0.0 .<= Pd[:,k] .<= (ps.shunt.P./ps.baseMVA)) # load served limits
        @constraint(m, [k=1:Ti], Pg[:,k] .<= ug[:,k].*Pg_max) # generator power limits upper
        @constraint(m, [k=1:Ti], ug[:,k].*Pg_min .<= Pg[:,k]) # generator power limits lower
        @constraint(m, [k=2:Ti], -ug[:,k].*ug[:,k-1].*RR .<= Pg[:,k]-Pg[:,k-1] .<= ug[:,k].*ug[:,k-1].*RR) # generator ramp rate
        @constraint(m, [l=-Tp+1:Ti], ug[:,k].*Pg_min .<= Pg[:,k]) # connecting ug to gon
        @constraint(m, [k=1:Ti], Ps_min .<= Ps[:,k] .<= Ps_max) # storage power flow
        @constraint(m, [j=2:Ti], E[:,j] .== (E[:,j-1] + ((dt/60) .*(Ps[:,j])))) # storage energy at next time step
        @constraint(m, [i=1:Ti], 0 .<= (E[:,i]) .<= E_max) # storage energy
        # set starting point (time at step 0 == k=1);
        @constraint(m, Pd[:,1] .== Pd1)
        @constraint(m, ug[:,1] .== ug1)
        @constraint(m, Pg[:,1] .== Pg1)
        @constraint(m, Ps[:,1] .== Ps1)
        @constraint(m, E[:,1] .== E1)
        # power balance
        @constraint(m, [k=1:Ti], 0 .== G_bus*Pg[:,k]+S_bus*Ps[:,k]-D_bus*Pd[:,k])
        # objective
        @objective(m, Max, 100*sum(sum(Pd)) + sum(sum(ug)));
        ## SOLVE! ##
        optimize!(m)
        sol_Pd=value.(Pd)[:,2]
        sol_Ps=value.(Ps)[:,2]
        sol_Pg=value.(Pg)[:,2]
        sol_E=value.(E)[:,2]
        sol_ug=value.(ug)[:,2]
        @assert abs(sum(ps.shunt.P.*ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg[gst]))<=2*tolerance
        @assert sum(ps.storage.E .< 0)==0
        dPd_star = (sol_Pd.*ps.baseMVA)./ps.shunt.P # % load served
        dPs_star = sol_Ps.*ps.baseMVA
        dPg_star = sol_Pg.*ps.baseMVA
        dE_star = sol_E.*ps.baseMVA
    end
    # add changes ps/psi structure
    ps.shunt.status = dPd_star;
    ps.storage.Ps = dPs_star;
    ps.storage.E = dE_star;
    ps.gen.Pg[gst] = dPg_star;
    return ps
end


function crisp_mh_rlopf!(ps,dt,t_win)
t = 0:dt:t_win;
Ti = size(t,1);
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
    RR = ps.gen[gst,:RampRateMWMin] .* dt ./ ps.baseMVA .* ps.gen[gst,:status]
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
    @variable(m, Pd[1:nd,1:Ti]) # demand
    @variable(m, Pg[1:ng,1:Ti]) # generation
    @variable(m, ug[1:ng,1:Ti], Bin) # generator j at time  k on == ug[j,k]=1, off == ug[j,k]=0
    @variable(m, Ps[1:ns,1:Ti]) # power flow into or out of storage (negative flow = charging)
    @variable(m, E[1:ns,1:Ti]) # energy level in battery
    @variable(m, Theta[1:n,1:Ti])
    # variable bounds constraints
    @constraint(m, [k=1:Ti], 0.0 .<= Pd[:,k] .<= ps.shunt.P./ps.baseMVA) # load served limits
    @constraint(m, [k=1:Ti], Pg[:,k] .<= ug[:,k].*Pg_max) # generator power limits upper
    @constraint(m, [k=1:Ti], ug[:,k].*Pg_min .<= Pg[:,k]) # generator power limits lower
    @constraint(m, [k=2:Ti], -ug[:,k].*ug[:,k-1].*RR .<= Pg[:,k]-Pg[:,k-1] .<= ug[:,k].*ug[:,k-1].*RR) # generator ramp rate lower
    @constraint(m, [k=1:Ti], Ps_min .<= Ps[:,k] .<= Ps_max) # storage power flow
    @constraint(m, [j=1:Ti-1], E[:,j+1] .== (E[:,j] + ((dt/60) .*Ps[:,j+1]))) # storage energy at next time step
    @constraint(m, [i=2:Ti], 0 .<= E[:,i] .<= E_max) # storage energy
    @constraint(m, [k=1:Ti], Theta[1,:] .== 0); # set first bus as reference bus: V angle to 0
    # power balance
    @constraint(m, [k=1:Ti], B*Theta[:,k] .== G_bus*Pg[:,k]+S_bus*Ps[:,k]-D_bus*Pd[:,k])
    # power flow limits
    @constraint(m, [k=1:Ti], -flow_max .<= Xinv.*(Theta[F,k] - Theta[T,k]) .<= flow_max)
    # objective
    @objective(m, Max, 100*sum(sum(Pd)) + sum(sum(ug)));
    ## SOLVE! ##
    optimize!(m)
    sol_Pd=value.(Pd)[:,2]
    sol_Ps=value.(Ps)[:,2]
    sol_Pg=value.(Pg)[:,2]
    sol_E=value.(E)[:,2]
    dE_star = sol_E.*ps.baseMVA
    dPd_star = sol_Pd.*ps.baseMVA./ps.shunt.P # % load served
    dPs_star = sol_Ps.*ps.baseMVA
    dPg_star = sol_Pg.*ps.baseMVA
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
    @variable(m, Pd[1:nd,1:Ti]) # demand
    @variable(m, Pg[1:ng,1:Ti]) # generation
    @variable(m, ug[1:ng,1:Ti], Bin) # generator j at time  k on == ug[j,k]=1, off == ug[j,k]=0
    @variable(m, Ps[1:ns,1:Ti]) # power flow into or out of storage (negative flow = charging)
    @variable(m, E[1:ns,1:Ti]) # energy level in battery
    # variable bounds constraints
    @constraint(m, [k=1:Ti], 0.0 .<= Pd[:,k] .<= ps.shunt.P./ps.baseMVA) # load served limits
    @constraint(m, [k=1:Ti], Pg[:,k] .<= ug[:,k].*Pg_max) # generator power limits upper
    @constraint(m, [k=1:Ti], ug[:,k].*Pg_min .<= Pg[:,k]) # generator power limits lower
    @constraint(m, [k=1:Ti], Ps_min .<= Ps[:,k] .<= Ps_max) # storage power flow
    @constraint(m, [j=1:Ti-1], E[:,j+1] .== (E[:,j] + ((dt/60) .*Ps[:,j+1]))) # storage energy at next time step
    @constraint(m, [i=2:Ti], 0 .<= E[:,i] .<= E_max) # storage energy
    # power balance
    @constraint(m, [k=1:Ti], 0 .== G_bus*Pg[:,k]+S_bus*Ps[:,k]-D_bus*Pd[:,k])
    # objective
    @objective(m, Max, 100*sum(sum(Pd)) + sum(sum(ug)));
    ## SOLVE! ##
    optimize!(m)
    sol_Pd=value.(Pd)[:,2]
    sol_Ps=value.(Ps)[:,2]
    sol_Pg=value.(Pg)[:,2]
    sol_E=value.(E)[:,2]
    dE_star = sol_E.*ps.baseMVA
    dPd_star = sol_Pd.*ps.baseMVA./ps.shunt.P # % load served
    dPs_star = sol_Ps.*ps.baseMVA
    dPg_star = sol_Pg.*ps.baseMVA
end
ps.shunt.status = dPd_star; #changes ps structure
ps.storage.Ps = dPs_star;
ps.storage.E = dE_star;
ps.gen.Pg[gst] = dPg_star; #changes ps structure
@assert abs(sum(ps.shunt.P.*ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg[gst]))<=2*tolerance
if !isempty(ps.storage.E)
    @assert sum(ps.storage.E .< -tolerance)==0
end
return ps
end
