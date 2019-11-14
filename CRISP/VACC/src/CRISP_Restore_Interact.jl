# moving horizon model of restoration process after blackout
# includes generator warm up and shut down times and storage
# still needs variable weather and load data to have appropriate hourly responses
using DataFrames;
using JuMP
using Gurobi
#using Cbc;
include("CRISP_network_gen.jl")
include("CRISP_interact.jl")
function crisp_Restoration_inter(ps,l_recovery_times,g_recovery_times,dt,t_window,t0,gen_on,comm,natg,nucp;load_cost=0)
    # constants
    tolerance = 10^(-6);
    if sum(load_cost)==0
        load_cost = ones(length(ps.shunt.P));
    end
    ti = t0;
    comm_count = zeros(length(ps.branch.t));
    comm_count[ps.branch.status .== 1] .= 100;
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
    EndTime = (t0+recTime+((maximum(ps.gen.minDownTimeHr)+maximum(ps.gen.minUpTimeHr))*60));
    Time = t0:dt:EndTime
    # find generator status
    ug = gen_on_off(ps,Time,t_window,gen_on,g_recovery_times)
    # find line status
    ul = line_stats(ps,Time,t_window,l_recovery_times)
    # varying load over the course of the optimization
    Pd_max = vary_load(ps,Time,t_window)
    # varying generation capacity over the optimization
    Pg_max = vary_gen_cap(ps,Time,t_window)
    # number of time steps within the time window
    t_win_step = Int64(t_window/dt);
    for i in 1:length(Time)
        # update time
        ti = Time[i]-t0;
        # remove failures as the recovery time is reached
        ps.branch.status[ti .>= l_recovery_times] .= 1;
        comm_count[ti .>= l_recovery_times] .= 100;
        ps.gen.status[ti .>= g_recovery_times] .= 1;
        # find the number of islands in ps
        subgraph = find_subgraphs(ps);# add Int64 here hide info here
        M = Int64(findmax(subgraph)[1]);
        ps_islands = build_islands(subgraph,ps)
        for j in 1:M
            psi = ps_subset(ps,ps_islands[j])
            i_subset = i:(i+t_win_step)
            ugi = ug[ps_islands[j].gen,i_subset]
            uli = ul[ps_islands[j].branch,i_subset]
            Pd_maxi = Pd_max[ps_islands[j].shunt,i_subset]
            Pg_maxi = Pg_max[ps_islands[j].gen,i_subset]
            crisp_mh_rlopf!(psi,dt,t_window,ugi,uli,Pd_maxi,Pg_maxi,load_cost[ps_islands[j].shunt])
            ps.gen.Pg[ps_islands[j].gen] = psi.gen.Pg
            ps.storage.Ps[ps_islands[j].storage] = psi.storage.Ps
            ps.storage.E[ps_islands[j].storage] = psi.storage.E
            ps.shunt.status[ps_islands[j].shunt] = psi.shunt.status
            #add_changes!(ps,psi,ps_islands[j]);
        end
        if comm
            if ti >= 4*60 #most communcation towers have batteries which have a capacity to cover from 4 to 24 hour
                for l in 1:length(ps.branch.f)
                    communication_interactions(ps,comm_count[l],ti)
                end
            end
        end
        if natg
            natural_gas_interactions(ps,)
        end
        if nucp
            nuclear_poissoning(ps)
        end
        # save current values
        cv.time .= ti;
        cv.load_shed .= sum(load_cost.*(Pd_max[:,i+1] - Pd_max[:,i+1].*ps.shunt.status));
        cv.perc_load_served .= (sum(load_cost.*Pd_max[:,i+1]) .- cv.load_shed)./sum(load_cost.*Pd_max[:,i+1]);
        cv.lines_out .= length(ps.branch.status) - sum(ps.branch.status);
        cv.gens_out .= length(ps.gen.status) - sum(ps.gen.status);
        append!(Restore,cv)
        @assert 10^(-4)>=abs(sum(Pd_max[:,i+1] .* ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg))
    end
    return Restore
end

function crisp_mh_rlopf!(ps,dt,t_win,ug,ul,Pd_max,Pg_max1,load_shed_cost)
    timeline = 0:dt:t_win
    Ti = size(timeline,1)
    # constants
    tolerance = 1e-4
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    bi = ps.bi
    # load data
    nd = size(ps.shunt,1)
    D = bi[ps.shunt.bus]
    D_bus = sparse(D,collect(1:nd),1.,n,nd)
    Pdmax = (Pd_max ./ ps.baseMVA)
    Pd1 = (Pd_max[:,1] ./ ps.baseMVA) .* ps.shunt.status
    if any(D.<1) || any(D.>n)
        error("Bad indices in shunt matrix")
    end
    # gen data
    ng = size(ps.gen.Pg,1)
    G = bi[ps.gen.bus]
    G_bus = sparse(G,collect(1:ng),1.,n,ng)
    Pg1 = (ps.gen.Pg ./ ps.baseMVA) .* ps.gen.status
    Pg_max = (Pg_max1 ./ ps.baseMVA)
    RR = ps.gen.RampRateMWMin .* dt ./ ps.baseMVA
    if any(G.<1) || any(G.>n)
        error("Bad indices in gen matrix")
    end
    # branch data
    nb = size(ps.branch,1)
    flow_max = ps.branch.rateA./ps.baseMVA # this could also be rateB
    # storage data
    min_bat_E = 0;
    ns = size(ps.storage,1)
    S = bi[ps.storage.bus]
    S_bus = sparse(S,collect(1:ns),1.,n,ns)
    Ps1 = (ps.storage.Ps ./ ps.baseMVA) .* ps.storage.status
    E1 = (ps.storage.E ./ ps.baseMVA) .* ps.storage.status
    E_max = (ps.storage.Emax ./ ps.baseMVA) .* ps.storage.status
    Ps_max = (ps.storage.Psmax ./ ps.baseMVA) .* ps.storage.status
    Ps_min = (ps.storage.Psmin ./ ps.baseMVA) .* ps.storage.status
    if any(S.<1) || any(S.>n)
        error("Bad indices in storage matrix")
    end
    # vector that depreciates the value of later elements in objective
    C_time = exp.(0:-1:(-Ti+1))';
    m = Model(with_optimizer(Gurobi.Optimizer))
    #m = Model(with_optimizer(Cbc.Optimizer))
    # variables
    @variable(m, Pd[1:nd, 1:Ti]) # demand
    @variable(m, Pg[1:ng, 1:Ti]) # generation
    @variable(m, Ps[1:ns, 1:Ti]) # power flow into or out of storage (negative flow = charging)
    @variable(m, E[1:ns, 1:Ti]) # energy level in battery
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
    @constraint(m, stPdcon[k=2:Ti], 0.0 .<= Pd[:,k] .<= Pdmax[:,k]) # load served limits
    @constraint(m, stPscon[k=2:Ti], Ps_min .<= Ps[:,k] .<= Ps_max) # storage power flow
    @constraint(m, stEPscon[k=2:Ti], E[:,k] .== (E[:,k-1] - ((dt/60) .* (Ps[:,k])))) # storage energy at next time step
    @constraint(m, stEcon[k=2:Ti], min_bat_E .<= (E[:,k]) .<= E_max) # storage energy
    @constraint(m, genPgucon[k=2:Ti], Pg[:,k] .<= ug[:,k] .* Pg_max[:,k]) # generator power limits upper
    @constraint(m, genPglcon[k=2:Ti], 0 .<= Pg[:,k]) # generator power limits lower
    if n > 1
        @variable(m, Theta[1:n, 1:Ti])
        @constraint(m, Theta[1,:] .== 0); # set first bus as reference bus: V angle to 0
        for k in 2:Ti
	  brst = falses(nb);
	  for b in 1:nb
	    if ul[b,k] == 1
              brst[b] = true;
	    end
	  end
            F = bi[ps.branch.f[brst]]
            T = bi[ps.branch.t[brst]]
            Xinv = (1 ./ ps.branch.X[brst])
            B = sparse(F,T,-Xinv,n,n) +
                sparse(T,F,-Xinv,n,n) +
                sparse(T,T,+Xinv,n,n) +
                sparse(F,F,+Xinv,n,n);
            #power balance
            @constraint(m, B*Theta[:,k] .== G_bus*Pg[:,k]+S_bus*Ps[:,k]-D_bus*Pd[:,k])
            # power flow limits
            @constraint(m, -flow_max[brst] .<= Xinv .* (Theta[F,k] - Theta[T,k]) .<= flow_max[brst])
        end
    else
        @constraint(m, PBnoTheta[k=2:Ti], 0.0 .== G_bus*Pg[:,k]+S_bus*Ps[:,k]-D_bus*Pd[:,k])
    end
    #ramping constraints
    @constraint(m, genPgRR[k=2:Ti], -RR .<= Pg[:,k-1] .- Pg[:,k] .<= RR) # generator ramp rate
    # objective
    @objective(m, Max, sum(Pd*C_time')) #TODO: constants
    ## SOLVE! ##
    optimize!(m)
    sol_Pd=value.(Pd)[:,2]
    sol_Ps=value.(Ps)[:,2]
    sol_Pg=value.(Pg)[:,2]
    sol_E=value.(E)[:,2]
    @assert abs(sum(Pd_max[:,2].*ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg))<=2*tolerance
    @assert sum(ps.storage.E .< -tolerance)==0
    dPd_star = (Vector(sol_Pd).*ps.baseMVA)./Pd_max[:,2] # % load served
    dPs_star = Vector(sol_Ps).*ps.baseMVA
    dPg_star = Vector(sol_Pg).*ps.baseMVA
    dE_star = Vector(sol_E).*ps.baseMVA
    # add changes ps/psi structure
    ps.shunt.status = dPd_star;
    ps.storage.Ps = dPs_star;
    ps.storage.E = dE_star;
    ps.gen.Pg = dPg_star;
    return ps
end

function gen_on_off(ps,Time,t_window,gen_on,gens_recovery_time)
    # constants
    tolerance = 1e-6
    ### collect the data that we will need ###
    dt = Time[2]-Time[1];
    Time = Time .- Time[1] .+ dt;
    ext_stps = Int64(t_window/dt);
    # gen data
    gs = (ps.gen.Pg .== 0);
    gst = (ps.gen.status .!= 1);
    g1 = (ps.gen.status .== 1);
    g2 = gen_on .& gs .& g1;
    ng = size(ps.gen.Pg,1)
    timedown = (ps.gen.minDownTimeHr .*60)
    timeup = (ps.gen.minDownTimeHr .*60)
    ug = falses(ng,length(Time)+ext_stps+1)
    gen_time = zeros(ng);
    gen_time[g2] .+= timedown[g2];
    gen_time[gs] .+= timeup[gs];
    gen_time[gst] .+= gens_recovery_time[gst];
    for t in 1:length(Time)+ext_stps+1
        gen_time .-= dt
        for g in 1:ng
            if gen_time[g] <= 0
                ug[g,t] = true;
            else
                ug[g,t] = false;
            end
        end
    end
    @assert sum(ug[:,end].!=true).==0
    return ug
end

function line_stats(ps,Time,t_window,l_recovery_time)
    ### collect the data that we will need ###
    dt = Time[2]-Time[1];
    Time = Time .- Time[1] .+ dt;
    ext_stps = Int64(t_window/dt);
    total_time_steps = length(Time)+ext_stps+1;
    # branch data
    bs = (ps.branch.status .== 0);
    nb = size(ps.branch,1);
    ul = falses(nb,length(Time)+ext_stps+1)
    branch_down_time = zeros(nb);
    branch_down_time[bs] .+= l_recovery_time[bs];
    for t in 1:total_time_steps
        branch_down_time .-= dt
        for b in 1:nb
            if branch_down_time[b] <= 0
                ul[b,t] = true;
            else
                ul[b,t] = false;
            end
        end
    end
    @assert sum(ul[:,end].!=true).==0
    return ul
end

function vary_load(ps,Time,t_window)
    Demand = CSV.File("data/solar+load/VT_hourly_load.csv") |> DataFrame
    DemandCap = maximum.(Demand.Load);
    TotD = Demand.Load ./ DemandCap .* sum(ps.shunt.P);
    PerD = ps.shunt.P ./ sum(ps.shunt.P);
    ### collect the data that we will need ###
    dt = Time[2]-Time[1];
    Time = Time .- Time[1] .+ dt;
    ext_stps = Int64(t_window/dt);
    # shunt data
    nd = size(ps.shunt.P,1);
    Pd = zeros(nd,length(Time)+ext_stps+1);
    if length(Time)+ext_stps+1 > length(TotD)
        n = Int64(round((length(Time)+ext_stps+1)/TotD))
        Pd = zeros(nd,n*length(TotD));
        for j in 0:n
            N = length(TotD)*j
            for t in 1:length(TotD)
                for d in 1:nd
                    Pd[d,N+t] = PerD[d]*TotD[N+t];
                end
            end
        end
    else
        for t in 1:length(Time)+ext_stps+1
            for d in 1:nd
                Pd[d,t] = PerD[d]*TotD[t];
            end
        end
    end
    @assert sum(Pd[:,:].<=0).==0
    return Pd
end

function vary_gen_cap(ps,Time,t_window)
    PercSolar1 = CSV.File("data/solar+load/NY_NE_1244665_solarPV_power_density.csv") |> DataFrame
    PowDen = PercSolar1.PowerDen;
    ### collect the data that we will need ###
    dt = Time[2]-Time[1];
    Time = Time .- Time[1] .+ dt;
    ext_stps = Int64(t_window/dt);
    # shunt data
    ng = size(ps.gen.Pg,1);
    Pg_max = zeros(ng,length(Time)+ext_stps+1);
    for t in 1:length(Time)+ext_stps+1
        for g in 1:ng
            if ps.gen.Fuel[g] .== "Solar"
                Pg_max[g,t] = ps.gen.Pmax[g].*PowDen[t];
            else
                Pg_max[g,t] = ps.gen.Pmax[g];
            end
        end
    end
    @assert sum(Pg_max[:,:].<0).==0
    return Pg_max
end
