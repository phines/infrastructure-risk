#set up packages
using CSV; using DataFrames; using SpecialFunctions;
include("CRISP_LSOPF_1.jl")
function RLSOPF!(totalp,ps,failures,recovery_times,Pd_max,load_cost=0)
    if load_cost==0
        load_cost = rand(length(ps.shunt[:P]));
    end
    rec_t = recovery_times[recovery_times.!=0];
    print(rec_t)
    times = sort(rec_t)
    load_shed = zeros(length(times)+1)
    # set load shed for the step just before restoration process
    load_shed[1] = totalp - sum(ps.shunt[:P]);
    for i = 1:length(times)
        T = times[i];
        # set failed branches to status 0
        failures[T.>=recovery_times] = 1;
        print(failures)
        # apply to network
        ps.branch[:,:status] = failures;
        # run the dcpf
        crisp_dcpf!(ps)
        # run lsopf
        (dPd, dPg) = crisp_rlopf(ps,Pd_max)
        # apply the results
        ps.gen[:Pg]  += dPg
        ps.shunt[:P] += dPd
        crisp_dcpf!(ps)
        # set load shed for this time step
        load_shed[i+1] = sum(Pd_max - ps.shunt[:P]);
    end
    times = [0;times];
    Restore = DataFrame(time = times, load_shed = load_shed)
    return Restore
end

function crisp_rlopf(ps,Pd_max)
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
    @constraint(m,-Pd.<=dPd.<=(Pd_max./ ps.baseMVA-Pd))
    @constraint(m,-Pg.<=dPg.<=(ps.gen[:Pmax]-Pg))
    @constraint(m,dTheta[1] == 0)
    # objective
    @objective(m,Max,sum(dPd) - 0.9*sum(dPg)) # serve as much load as possible
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
