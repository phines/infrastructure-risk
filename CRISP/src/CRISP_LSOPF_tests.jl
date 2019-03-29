using JuMP
using Clp
using SparseArrays
using LinearAlgebra

#export run_dcpf

function crisp_dcpf!(ps)
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
    # gen data
    ng = size(ps.gen,1)
    G = bi[ps.gen[:bus]]
    Pg = ps.gen[:Pg] ./ ps.baseMVA .* ps.gen[:status]
    # branch data
    brst = (ps.branch[:status].==1)
    F = bi[ps.branch[brst,:f]]
    T = bi[ps.branch[brst,:t]]
    Xinv = (1 ./ ps.branch[brst,:X])
    Bdc = sparse(F,T,-Xinv,n,n) + sparse(T,F,-Xinv,n,n) +
          sparse(T,T,+Xinv,n,n) + sparse(F,F,+Xinv,n,n)
    # bus injection
    Pbus = Array(sparse(G,fill(1,ng),Pg,n,1) - sparse(D,fill(1,nd),Pd,n,1))
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

function crisp_lsopf!(ps)
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    if n>1
        bi = sparse(ps.bus.id,fill(1,n),collect(1:n)) # helps us to find things
        # load data
        nd = size(ps.shunt,1)
        D = bi[ps.shunt.bus]
        Pd = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
        # gen data
        ng = size(ps.gen,1)
        G = bi[ps.gen.bus]
        Pg = ps.gen.Pg ./ ps.baseMVA .* ps.gen.status
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
        @variable(m,dTheta[1:n])
        # variable bounds
        @constraint(m,-Pd.<=dPd.<=0)
        @constraint(m,-Pg.<=dPg.<=0)
        @constraint(m,dTheta[1] == 0)
        # objective
        @objective(m,Max,sum(dPd)) # serve as much load as possible
        # mapping matrix to map loads/gens to buses
        M_D = sparse(D,1:nd,1.0,n,nd)
        M_G = sparse(G,1:ng,1.0,n,ng)
        # Power balance equality constraint
        @constraint(m,B*dTheta .== M_G*dPg - M_D*dPd)
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
    else
        if (!isempty(ps.gen) && isempty(ps.shunt)) || (isempty(ps.gen) && !isempty(ps.shunt))
            deltaPg = -(ps.gen.Pg ./ ps.baseMVA .* ps.gen.status);
            deltaPd = -(ps.shunt.P ./ ps.baseMVA .* ps.shunt.status);
            deltaPg_star = deltaPg.*ps.baseMVA;
            deltaPd_star = deltaPd.*ps.baseMVA;
            ps.gen.Pg  += deltaPg_star;
            ps.shunt.P += deltaPd_star;
        elseif !isempty(ps.gen) && !isempty(ps.shunt)
            Pd = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
            Pg = ps.gen.Pg ./ ps.baseMVA .* ps.gen.status
            Pg_cap = ps.gen.Pmax ./ ps.baseMVA .* ps.gen.status
            if sum(Pg_cap) >= sum(Pd)
                deltaPd = 0.0;
                deltaPg = sum(Pd)-sum(Pg);
            else
                deltaPd = sum(Pd)-sum(Pg_cap);
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
            ps.shunt.P = deltaPd_star;
            ps.gen.Pg  = deltaPg_star;
        else
            ps.shunt.P = ps.shunt.P.*0.0;
            ps.gen.Pg  = ps.gen.Pg.*0.0;
        end
    end
    return ps
end
#end
