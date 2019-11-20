using JuMP
using Clp
using Cbc
using Gurobi
using SparseArrays
using LinearAlgebra

function crisp_dcpf_g_s!(ps)
    # constants
    tolerance = 1e-6
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    bi = sparse(ps.bus.id,fill(1,n),collect(1:n)) # helps us to find things
    # load data
    nd = size(ps.shunt,1)
    D = bi[ps.shunt.bus]
    Pd = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
    if any(D.<1) || any(D.>n)
        error("Bad indices in shunt matrix")
    end
    Pd_bus = Array(sparse(D,ones(size(D)),Pd,n,1))
    # gen data
    gst = (ps.gen.status.==1)
    ng = size(ps.gen[gst,:],1)
    G = bi[ps.gen.bus[gst]]
    Pg = ps.gen.Pg[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    if any(G.<1) || any(G.>n)
        error("Bad indices in gen matrix")
    end
    Pg_bus = Array(sparse(G,ones(size(G)),Pg,n,1))
    # storage data
    S = bi[ps.storage.bus];
    Ps = ps.storage.Ps ./ ps.baseMVA .* ps.storage.status
    if any(S.<1) || any(S.>n)
        error("Bad indices in sstorage matrix")
    end
        Ps_bus = Array(sparse(S,ones(size(S)),Ps,n,1))
    # branch data
    brst = (ps.branch.status .== 1)
    F = bi[ps.branch.f[brst]]
    T = bi[ps.branch.t[brst]]
    Xinv = (1 ./ ps.branch.X[brst])
    Bdc = sparse(F,T,-Xinv,n,n) + sparse(T,F,-Xinv,n,n) +
          sparse(T,T,+Xinv,n,n) + sparse(F,F,+Xinv,n,n)
    #find reference bus
    if isempty(ps.bus.id[ps.bus.bus_type.==3]) # note in Pavan's ps structure I beleive it's called 'kind' not 'bus_type'
        if (isempty(ps.gen.Pg[gst]))  || (isempty(ps.shunt) && isempty(ps.storage))
            ps.branch.Pf[brst] .= 0
            ps.branch.Pt[brst] .= 0
            ps.branch.Qf[brst] .= 0
            ps.branch.Qt[brst] .= 0
            if !isempty(ps.gen.Pg[gst])
                ps.gen.Pg[gst] .= 0.0
            end
            if !isempty(ps.shunt)
                ps.shunt.P .= 0.0
            end
            if !isempty(ps.storage)
                ps.storage.Ps .= 0.0
            end
            return ps #TODO
        elseif !isempty(ps.gen.Pg[gst])
            maxGen = findmax(ps.gen.Pg)[2]
            busID = ps.gen[maxGen,:bus];
            isref = (busID.==ps.bus.id)
            nonref = .~isref
        end
    elseif sum(ps.gen.status[ps.gen.bus.==ps.bus.id[ps.bus.bus_type.==3]])==0
        if isempty(ps.gen[gst]) || (isempty(ps.shunt) && isempty(ps.storage))
            ps.branch.Pf[brst] .= 0
            ps.branch.Pt[brst] .= 0
            ps.branch.Qf[brst] .= 0
            ps.branch.Qt[brst] .= 0
            if !isempty(ps.gen[gst])
                ps.gen.Pg[gst] .= 0.0
            end
            if !isempty(ps.shunt)
                ps.shunt.P .= 0.0
            end
            if !isempty(ps.storage)
                ps.storage.Ps .= 0.0
            end
            return ps #TODO
        else
            maxGen = findmax(ps.gen.Pg)[2]
            busID = ps.gen[maxGen,:bus];
            isref = (busID.==ps.bus.id)
            nonref = .~isref
        end
    else
        if isempty(ps.gen.Pg[gst]) || (isempty(ps.shunt) && isempty(ps.storage))
            ps.branch.Pf[brst] .= 0
            ps.branch.Pt[brst] .= 0
            ps.branch.Qf[brst] .= 0
            ps.branch.Qt[brst] .= 0
            if !isempty(ps.gen.Pg[gst])
                ps.gen.Pg[gst] .= 0.0
            end
            if !isempty(ps.shunt)
                ps.shunt.P .= 0.0
            end
            if !isempty(ps.storage)
                ps.storage.Ps .= 0.0
            end
            return ps #TODO
        else
            isref = (ps.bus.bus_type.==3)
            nonref = .~isref
        end
    end
    # bus injection
    Pbus = Pg_bus+Ps_bus-Pd_bus;#Array(sparse(G,fill(1,ng),Pg,n,1) - sparse(D,fill(1,nd),Pd,n,1))
    # angles
    theta = zeros(n)
    Bsub = Bdc[nonref,nonref]
    Psub = Pbus[nonref]
    tsub = Bsub\Psub
    theta[nonref] = tsub
    # record the results to the bus matrix
    ps.bus.Va = theta .* (180.0 / pi)
    ps.bus.Vm = ones(n)
    # compute/record the power flows
    Pf_pu = Xinv .* (theta[F] - theta[T])
    ps.branch.Pf[brst] = +Pf_pu.*ps.baseMVA
    ps.branch.Pt[brst] = -Pf_pu.*ps.baseMVA
    ps.branch.Qf[brst] .= 0
    ps.branch.Qt[brst] .= 0
    # fix the generation at the slack bus
    mismatch = sum(Pbus)
    if abs(mismatch)>tolerance
        refbusid = ps.bus.id[isref]
        is_refgen = (ps.gen.bus .==refbusid)
        if sum(ps.gen.status[is_refgen].!=1) >=1
            is_refgen[.!gst] .= false;
        end
        if sum(is_refgen) >= 1
            println("Multiple Gen on ref bus, splitting mismatch among them.")
            ps.gen.Pg[is_refgen] .-= (mismatch.*ps.baseMVA)/sum(is_refgen)
        elseif sum(is_refgen) == 1
            ps.gen.Pg[is_refgen] .-= (mismatch.*ps.baseMVA)
        else
            println(ps.gen)
            println(is_refgen)
            error("No reference generator")
        end
    end
    # check the mismatch
    mis_check = sum(ps.gen.Pg[gst].*ps.gen.status[gst]) + sum(ps.storage.Ps.*ps.storage.status) - sum(ps.shunt.P.*ps.shunt.status)
    if abs(mis_check)>tolerance
        println(mismatch)
        println(Pbus)
        print("Mismatch = ")
        println(mis_check)
        error("Mismach error in crisp_dcpf")
    end
    # return the resulting system
    return ps
end

function crisp_lsopf_g1_s!(ps,dt)
    # constants
    tolerance = 1e-4
    load_shed_cost_C = 100
    min_bat_E = 0.01;
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    bi = sparse(ps.bus.id,fill(1,n),collect(1:n)) # helps us to find things
    # load data
    nd = size(ps.shunt,1)
    load_shed_cost = load_shed_cost_C.*ones(nd)
    D = bi[ps.shunt.bus]
    D_bus = sparse(D,collect(1:nd),1.,n,nd);
    Pd = ps.shunt.P  ./ ps.baseMVA
    if any(D.<1) || any(D.>n)
        error("Bad indices in shunt matrix")
    end
    #Pd_bus = Array(sparse(D,ones(size(D)),Pd,n,1))
    # gen data
    gst = (ps.gen.status.==1)
    ng = length(ps.gen.Pg[gst])
    G = bi[ps.gen.bus[gst]]
    G_bus = sparse(G,collect(1:ng),1.,n,ng);
    Pg = ps.gen.Pg[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    RR = ps.gen.RampRateMWMin[gst] .* dt ./ ps.baseMVA
    Pg_max = ps.gen.Pmax[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    Pg_min = ps.gen.Pmin[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    if any(G.<1) || any(G.>n)
        error("Bad indices in gen matrix")
    end
    #Pg_bus = Array(sparse(G,ones(size(G)),Pg,n,1))
    # storage data
    ns = size(ps.storage,1)
    S = bi[ps.storage.bus];
    S_bus = sparse(S,collect(1:ns),1.,n,ns);
    Ps = ps.storage.Ps ./ ps.baseMVA .* ps.storage.status
    E = ps.storage.E ./ ps.baseMVA .* ps.storage.status
    E_max = ps.storage.Emax ./ ps.baseMVA .* ps.storage.status
    Ps_max = ps.storage.Psmax ./ ps.baseMVA .* ps.storage.status
    Ps_min = ps.storage.Psmin ./ ps.baseMVA .* ps.storage.status
    if any(S.<1) || any(S.>n)
        error("Bad indices in storage matrix")
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
    m = Model(with_optimizer(Gurobi.Optimizer))
    #m = Model(with_optimizer(Cbc.Optimizer))
        # variables
        @variable(m,dPd[1:nd])
        @variable(m,dPs[1:ns])
        @variable(m,ndPg[1:ng])
        @variable(m,pdPg[1:ng])
        @variable(m,ug[1:ng],Bin)
        # variable bounds
        @constraint(m,-Pd .<= dPd .<= 0)
        @constraint(m, Ps_min .<= Ps+dPs .<= Ps_max)
        @constraint(m, 0 .<= E - (Ps+dPs).*dt/60 .<= E_max)
        @constraint(m, ug.*(Pg_min) .<= Pg+ndPg+pdPg)
        @constraint(m, ug.*(Pg_max) .>= Pg+pdPg+ndPg)
        for g in 1:ng
            if ug[g] .== 1
                #ramp rate constraints
                @constraint(m, -RR[g] .<= ndPg[g]+pdPg[g] .<= RR[g])
            else
            end
        end
        @constraint(m, pdPg .>= 0)
        @constraint(m, ndPg .<= 0)
        if n>1
            @variable(m,dTheta[1:n])
            @constraint(m,dTheta[1] == 0)
            # Power balance equality constraint
            @constraint(m,B*dTheta .== G_bus*ndPg+G_bus*pdPg+S_bus*dPs-D_bus*dPd);#M_G*dPg - M_D*dPd)
            # Power flow constraints
            @constraint(m,-flow_max .<= flow0 + Xinv.*(dTheta[F] - dTheta[T]) .<= flow_max)
        else
            # Power balance equality constraint
            @constraint(m,sum(Pg+ndPg+pdPg)+sum(Ps+dPs)-sum(Pd+dPd) == 0);
        end
        # objective
        @objective(m,Max,sum(load_shed_cost.*dPd)+sum(ug)) # serve as much load as possible and keep as much generation as possible
    ### solve the model ###
    optimize!(m);
    # collect/return the outputs
    sol_dPd=value.(dPd)
    sol_dPs=value.(dPs)
    sol_ndPg=value.(ndPg)
    sol_pdPg=value.(pdPg)
    dPd_star = sol_dPd.*ps.baseMVA./ps.shunt.P
    dPs_star = sol_dPs.*ps.baseMVA
    dPg_star = (sol_pdPg+sol_ndPg).*ps.baseMVA
    ps.shunt.status += dPd_star; #changes ps structure
    ps.storage.Ps += dPs_star;
    ps.storage.E += ps.storage.Ps.*dt;
    ps.gen.Pg[gst] += dPg_star; #changes ps structure
    #adding criteria that should produce errors if incorrect.
    @assert abs(sum(ps.shunt.P.*ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg[gst]))<=2*tolerance
    #@assert sum(ps.storage.E .< 0)==0
    #@assert 0.<=ps.shunt.P
    return ps
end


## CRISP_LSOPF_gen_stor.jl
function crisp_lsopf_g_s!(ps,dt)
t = [0; dt];
Ti = size(t,1);
# constants
tolerance = 1e-5
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
    g1 = (ps.gen.status .== 1);
    g2 = (ps.gen.Pg .!=0);
    g3 = (ps.gen.time_on .!= 0)
    gst = g1 .| g2;
    gf = .!gst .& g3
    ps.gen.time_off[gf] .+= dt/60;
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
    @objective(m, Max, 100*sum(sum(Pd)) + sum(sum(ug)));
    ## SOLVE! ##
    optimize!(m)
    sol_Pd=value.(Pd)
    sol_Ps=value.(Ps)
    sol_Pg=value.(Pg)
    sol_E=value.(E)
    sol_ug=value.(ug)
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
    g1 = (ps.gen.status .== 1);
    g2 = (ps.gen.Pg .!=0);
    g3 = (ps.gen.time_on .!= 0)
    gst = g1 .| g2;
    gf = gst .& g3
    ps.gen.time_off[gf] .+= dt/60;
    ng = size(ps.gen[gst,:Pg],1)
    G = bi[ps.gen[gst,:bus]]
    G_bus = sparse(G,collect(1:ng),1.,n,ng);
    Pg1 = ps.gen[gst,:Pg] ./ ps.baseMVA .* ps.gen[gst,:status]
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
    dE_star = sol_E.*ps.baseMVA
    dPd_star = sol_Pd.*ps.baseMVA./ps.shunt.P # % load served
    dPs_star = sol_Ps.*ps.baseMVA
    dPg_star = sol_Pg.*ps.baseMVA
end
ps.shunt.status = dPd_star; #changes ps structure
ps.storage.Ps = dPs_star;
ps.storage.E = dE_star;
ps.gen.Pg[gst] = dPg_star; #changes ps structure
#ps.gen.Pg[gst][sol_ug.==0] .= 0.0;
ps.gen.time_off[gst][sol_ug.==0] .+= dt/60;
@assert abs(sum(ps.shunt.P.*ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg[gst]))<=2*tolerance
if !isempty(ps.storage.E)
    @assert sum(ps.storage.E .< -tolerance)==0
end
return ps
end

function crisp_opf_initiate!(ps,dt)
t = [0; dt];
Ti = size(t,1);
# constants
tolerance = 1e-5
### collect the data that we will need ###
# bus data
n = size(ps.bus,1) # the number of buses
if n>1
    bi = ps.bi;
    # load data
    nd = size(ps.shunt,1)
    D = bi[ps.shunt[:bus]]
    D_bus = sparse(D,collect(1:nd),1.,n,nd);
    Pd1 = ps.shunt[:P] ./ ps.baseMVA;
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
    @variable(m, Pg[1:ng]) # generation
    @variable(m, ug[1:ng], Bin) # generator j at time  k on == ug[j,k]=1, off == ug[j,k]=0
    @variable(m, Ps[1:ns]) # power flow into or out of storage (negative flow = charging)
    @variable(m, E[1:ns]) # energy level in battery
    @variable(m, Theta[1:n])
    # variable bounds constraints
    @constraint(m, Pg .<= ug.*Pg_max) # generator power limits upper
    @constraint(m, ug.*Pg_min .<= Pg) # generator power limits lower
    @constraint(m, Ps_min .<= Ps .<= Ps_max) # storage power flow
    @constraint(m, E .== (E1 + (dt/60) .*(Ps))) # storage energy at next time step
    @constraint(m, 0 .<= E .<= E_max) # storage energy
    @constraint(m, Theta[1] .== 0); # set first bus as reference bus: V angle to 0
    # power balance
    @constraint(m, B*Theta .== G_bus*Pg+S_bus*Ps-D_bus*Pd1)
    # power flow limits
    @constraint(m, -flow_max .<= Xinv.*(Theta[F] - Theta[T]) .<= flow_max)
    # objective
    @objective(m, Max, -sum(sum(Ps)) + sum(sum(ug)));
    ## SOLVE! ##
    optimize!(m)
    sol_Ps=value.(Ps)
    sol_Pg=value.(Pg)
    sol_E=value.(E)
    sol_ug=value.(ug)
    dE_star = sol_E.*ps.baseMVA
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
    g1 = (ps.gen.status .== 1);
    g2 = (ps.gen.Pg .!=0);
    g3 = (ps.gen.time_on .== 0)
    ps.gen.time_off[!g2 && g3] += dt./60;
    gst = g1 && g2;
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
    @variable(m, Pg[1:ng]) # generation
    @variable(m, ug[1:ng], Bin) # generator j at time  k on == ug[j,k]=1, off == ug[j,k]=0
    @variable(m, Ps[1:ns]) # power flow into or out of storage (negative flow = charging)
    @variable(m, E[1:ns]) # energy level in battery
    # variable bounds constraints
    @constraint(m, Pg .<= ug.*Pg_max) # generator power limits upper
    @constraint(m, ug.*Pg_min .<= Pg) # generator power limits lower
    @constraint(m, Ps_min .<= Ps .<= Ps_max) # storage power flow
    @constraint(m, E .== (E1 + (dt/60) .*(Ps))) # storage energy at next time step
    @constraint(m, 0 .<= E .<= E_max) # storage energy
    # power balance
    @constraint(m, 0 .== G_bus*Pg+S_bus*Ps-D_bus*Pd1)
    # objective
    @objective(m, Max, 100*sum(sum(Pd)) + sum(sum(ug)));
    ## SOLVE! ##
    optimize!(m)
    sol_Ps=value.(Ps)
    sol_Pg=value.(Pg)
    sol_E=value.(E)
    sol_ug=value.(ug)
    dE_star = sol_E.*ps.baseMVA
    dPs_star = sol_Ps.*ps.baseMVA
    dPg_star = sol_Pg.*ps.baseMVA
end
ps.storage.Ps = dPs_star;
ps.storage.E = dE_star;
ps.gen.Pg[gst] = dPg_star; #changes ps structure
ps.gen.Pg[sol_ug.==0] .= 0.0;
ps.gen.time_off[sol_ug.==0] .+= ps.gen.minDownTimeHr[sol_ug.==0];
@assert abs(sum(ps.shunt.P.*ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg[gst]))<=2*tolerance
return ps
end
