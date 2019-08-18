using JuMP
using Clp
using Cbc
#using Gurobi
using SparseArrays
using LinearAlgebra

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
