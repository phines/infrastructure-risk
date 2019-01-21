#module CRISP_LSOPF
using JuMP
using Clp
#using SparseArrays
#using LinearAlgebra

#export run_dcpf
function crisp_dcpf!(ps)
    # constants
    eps = 1e-6
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    bi = sparse(ps.bus[:id],fill(1,n),collect(1:n)) # helps us to find things
    isref = (ps.bus[:bus_type].==3)
    nonref = .~isref
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
    tsub = Bsub\Psub # get an error when network is split into pieces
    theta[nonref] = tsub
    # record the results to the bus matrix
    ps.bus[9] = theta .* (180.0 / pi)
    ps.bus[8] = ones(n)
    # compute/record the power flows
    Pf_pu = Xinv .* (theta[F] - theta[T])
    ps.branch[brst,:Pf] = +Pf_pu.*ps.baseMVA
    ps.branch[brst,:Pt] = -Pf_pu.*ps.baseMVA
    ps.branch[brst,:Qf] .= 0
    ps.branch[brst,:Qt] .= 0
    # fix the generation at the slack bus
    mismatch = sum(Pbus)
    if abs(mismatch)>eps
        refbusid = ps.bus[isref,:id]
        is_refgen = (ps.gen[:bus].==refbusid)
        if sum(is_refgen) != 1
            error("Must be exactly one ref generator")
        end
        ps.gen[is_refgen,:Pg] -= (mismatch.*ps.baseMVA)
    end
    # check the mismatch
    mis_check = sum(ps.gen[:,:Pg].*ps.gen[:,:status]) - sum(ps.shunt[:,:P].*ps.shunt[:,:status]);
    if abs(mis_check)>eps
        println(mismatch)
        println(Pbus)
        print("Mismatch = ")
        println(mis_check)
        error("Mismach error in crisp_dcpf")
    end
    # return the resulting system
    return ps
end

function crisp_lsopf(ps)
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    bi = sparse(ps.bus[:id],fill(1,n),collect(1:n)) # helps us to find things
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
    flow0 = ps.branch[brst,:Pf]./ps.baseMVA
    flow_max = ps.branch[brst,:rateA]./ps.baseMVA # this could also be rateB
    Xinv = (1 ./ ps.branch[brst,:X])
    B = sparse(F,T,-Xinv,n,n) +
        sparse(T,F,-Xinv,n,n) +
        sparse(T,T,+Xinv,n,n) +
        sparse(F,F,+Xinv,n,n)
    ### Build the optimization model ###
    m = Model(solver = ClpSolver())
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
    solve(m)
    # collect/return the outputs
    dPd_star = getvalue(dPd).*ps.baseMVA
    dPg_star = getvalue(dPg).*ps.baseMVA
    return (dPd_star, dPg_star)
end

#end
