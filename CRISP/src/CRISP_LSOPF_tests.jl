using JuMP
using Clp
using SparseArrays
using LinearAlgebra

#export run_dcpf
function crisp_dcpf1!(ps)
    # constants
    tolerance = 1e-6
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    bi = sparse(ps.bus.id,fill(1,n),collect(1:n)) # helps us to find things
    if isempty(ps.bus.id[ps.bus.bus_type.==3]) # note in Pavan's ps structure I beleive it's called 'kind' not 'bus_type'
        if isempty(ps.gen) || isempty(ps.shunt)
            return ps #TODO: is this working?
        else
            maxGen = findmax(ps.gen.Pg)[2]
            busID = ps.gen[maxGen,:bus];
            isref = (busID.==ps.bus.id)
            nonref = .~isref
        end
    else
        isref = (ps.bus.bus_type.==3)
        nonref = .~isref
    end
    #isref = (ps.bus[:bus_type].==3)
    #nonref = .~isref
    # load data
    nd = size(ps.shunt,1)
    D = bi[ps.shunt[:bus]]
    Pd = ps.shunt[:P] ./ ps.baseMVA .* ps.shunt[:status]
    if any(D.<1) || any(D.>n)
        error("Bad indices in shunt matrix")
    end
    Pd_bus = Array(sparse(D,ones(size(D)),Pd,n,1))
    # gen data
    ng = size(ps.gen,1)
    G = bi[ps.gen[:bus]]
    Pg = ps.gen[:Pg] ./ ps.baseMVA .* ps.gen[:status]
    if any(G.<1) || any(G.>n)
        error("Bad indices in gen matrix")
    end
    Pg_bus = Array(sparse(G,ones(size(G)),Pg,n,1))
    # branch data
    brst = (ps.branch[:status].==1)
    F = bi[ps.branch[brst,:f]]
    T = bi[ps.branch[brst,:t]]
    Xinv = (1 ./ ps.branch[brst,:X])
    Bdc = sparse(F,T,-Xinv,n,n) + sparse(T,F,-Xinv,n,n) +
          sparse(T,T,+Xinv,n,n) + sparse(F,F,+Xinv,n,n)
    # bus injection
    Pbus = Pg_bus-Pd_bus;#Array(sparse(G,fill(1,ng),Pg,n,1) - sparse(D,fill(1,nd),Pd,n,1))
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
        refbusid = ps.bus[isref,:id]
        is_refgen = (ps.gen[:bus].==refbusid)
        if sum(is_refgen) != 1
            error("Must be exactly one ref generator")
        end
        ps.gen.Pg[is_refgen] .-= (mismatch.*ps.baseMVA)
    end
    # check the mismatch
    mis_check = sum(ps.gen.Pg.*ps.gen.status) - sum(ps.shunt.P.*ps.shunt.status)
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

function crisp_lsopf1!(ps)
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
        Pd_bus = Array(sparse(D,ones(size(D)),Pd,n,1))
        # gen data
        ng = size(ps.gen,1)
        G = bi[ps.gen[:bus]]
        G_bus = sparse(G,collect(1:ng),1.,n,ng);
        Pg = ps.gen[:Pg] ./ ps.baseMVA .* ps.gen[:status]
        if any(G.<1) || any(G.>n)
            error("Bad indices in gen matrix")
        end
        Pg_bus = Array(sparse(G,ones(size(G)),Pg,n,1))
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
        m = Model(with_optimizer(Clp.Optimizer))
        # variables
        @variable(m,dPd[1:nd])
        @variable(m,dPg[1:ng])
        @variable(m,dPd_bus[1:n])
        @variable(m,dPg_bus[1:n])
        @variable(m,dTheta[1:n])
        # variable bounds
        @constraint(m,-Pd.<=dPd.<=0)
        @constraint(m,-Pg.<=dPg.<=0)
        @constraint(m,-Pd_bus.<=dPd_bus.<=0)
        @constraint(m,-Pg_bus.<=dPg_bus.<=0)
        @constraint(m,dTheta[1] == 0)
        #coupling constraint
        @constraint(m,dPd_bus .== D_bus*dPd) #TODO only works when there are not multiple loads per bus...need ==Array(sparse(D,ones(size(D)),dPd,n,1)))
        @constraint(m,dPg_bus .== G_bus*dPg) #TODO only works when there are not multiple gens per bus...need ==Array(sparse(G,ones(size(G)),dPg,n,1)))
        # objective
        @objective(m,Max,sum(dPd)) # serve as much load as possible
        # mapping matrix to map loads/gens to buses
        #M_D = (sparse(D,ones(size(D)),1.0,n,1));#sparse(D,1:nd,1.0,n,nd)
        #M_G = (sparse(G,ones(size(G)),1.0,n,1));#sparse(G,1:ng,1.0,n,ng)
        # Power balance equality constraint
        @constraint(m,B*dTheta .== dPg_bus-dPd_bus);#M_G*dPg - M_D*dPd)
        # Power flow constraints
        @constraint(m,-flow_max .<= flow0 + Xinv.*(dTheta[F] - dTheta[T]) .<= flow_max)
        ### solve the model ###
        optimize!(m);
        # collect/return the outputs
        sol_dPd=value.(dPd)
        sol_dPg=value.(dPg)
        dPd_star = sol_dPd.*ps.baseMVA
        dPg_star = sol_dPg.*ps.baseMVA
        ps.shunt.P += dPd_star; #changes ps structure
        ps.gen.Pg += dPg_star; #changes ps structure

    return ps
end
#end