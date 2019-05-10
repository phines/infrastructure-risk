using JuMP
using Clp
using Cbc
using SparseArrays
using LinearAlgebra

#export run_dcpf
function crisp_dcpf_g!(ps)
    # constants
    tolerance = 1e-6
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    bi = sparse(ps.bus.id,fill(1,n),collect(1:n)) # helps us to find things
    # load data
    nd = size(ps.shunt,1)
    D = bi[ps.shunt[:bus]]
    Pd = ps.shunt[:P] ./ ps.baseMVA .* ps.shunt[:status]
    if any(D.<1) || any(D.>n)
        error("Bad indices in shunt matrix")
    end
    Pd_bus = Array(sparse(D,ones(size(D)),Pd,n,1))
    # gen data
    gst = (ps.gen[:status].==1)
    ng = size(ps.gen[gst,:],1)
    G = bi[ps.gen[gst,:bus]]
    Pg = ps.gen[gst,:Pg] ./ ps.baseMVA .* ps.gen[gst,:status]
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
    #find reference bus
    if isempty(ps.bus.id[ps.bus.bus_type.==3]) # note in Pavan's ps structure I beleive it's called 'kind' not 'bus_type'
        if isempty(ps.gen[gst,:]) || isempty(ps.shunt)
            ps.branch.Pf[brst] .= 0
            ps.branch.Pt[brst] .= 0
            ps.branch.Qf[brst] .= 0
            ps.branch.Qt[brst] .= 0
            if isempty(ps.shunt) && !isempty(ps.gen[gst,:])
                ps.gen.Pg[gst] .= 0.0
            elseif !isempty(ps.shunt) && isempty(ps.gen[gst,:])
                ps.shunt.P .= 0.0
            end
            return ps #TODO
        else
            maxGen = findmax(ps.gen.Pg[gst])[2]
            busID = ps.gen[maxGen,:bus];
            isref = (busID.==ps.bus.id)
            nonref = .~isref
        end
    elseif !isempty(ps.bus.id[ps.bus.bus_type.==3]) && ps.gen.status[ps.gen.bus.==ps.bus.id[ps.bus.bus_type.==3]]!=1
        if isempty(ps.gen[gst,:]) || isempty(ps.shunt)
            ps.branch.Pf[brst] .= 0
            ps.branch.Pt[brst] .= 0
            ps.branch.Qf[brst] .= 0
            ps.branch.Qt[brst] .= 0
            if isempty(ps.shunt) && !isempty(ps.gen[gst,:])
                ps.gen.Pg[gst] .= 0.0
            elseif !isempty(ps.shunt) && isempty(ps.gen[gst,:])
                ps.shunt.P .= 0.0
            end
            return ps #TODO
        else
            maxGen = findmax(ps.gen.Pg[gst])[2]
            busID = ps.gen[maxGen,:bus];
            isref = (busID.==ps.bus.id)
            nonref = .~isref
        end
    else
        isref = (ps.bus.bus_type.==3)
        nonref = .~isref
    end
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
        if sum(ps.gen.status[is_refgen].!=1) >=1
            is_refgen[.!gst] .= false;
        end
        if sum(is_refgen) >= 1
            println("Multiple Gen on ref bus, splitting mismatch among them.")
            ps.gen.Pg[is_refgen] .-= (mismatch.*ps.baseMVA)/sum(is_refgen)
        elseif sum(is_refgen) == 1
            ps.gen.Pg[is_refgen] .-= (mismatch.*ps.baseMVA)
        else
            println(is_refgen)
            error("No reference generator")
        end
    end
    # check the mismatch
    mis_check = sum(ps.gen.Pg[gst].*ps.gen.status[gst]) - sum(ps.shunt.P.*ps.shunt.status)
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

function crisp_lsopf_g!(ps)
    # constants
    tolerance = 1e-6
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    if n>1
        bi = sparse(ps.bus.id,fill(1,n),collect(1:n)) # helps us to find things
        # load data
        nd = size(ps.shunt,1)
        D = bi[ps.shunt[:bus]]
        D_bus = sparse(D,collect(1:nd),1.,n,nd);
        Pd = ps.shunt[:P] ./ ps.baseMVA .* ps.shunt[:status]
        if any(D.<1) || any(D.>n)
            error("Bad indices in shunt matrix")
        end
        #Pd_bus = Array(sparse(D,ones(size(D)),Pd,n,1))
        # gen data
        gst = (ps.gen[:status].==1)
        ng = size(ps.gen[gst,:],1)
        G = bi[ps.gen[gst,:bus]]
        G_bus = sparse(G,collect(1:ng),1.,n,ng);
        Pg = ps.gen[gst,:Pg] ./ ps.baseMVA .* ps.gen[gst,:status]
        Pg_max = ps.gen[gst,:Pmax] ./ ps.baseMVA .* ps.gen[gst,:status]
        Pg_min = ps.gen[gst,:Pmin] ./ ps.baseMVA .* ps.gen[gst,:status]
        if any(G.<1) || any(G.>n)
            error("Bad indices in gen matrix")
        end
        #Pg_bus = Array(sparse(G,ones(size(G)),Pg,n,1))
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
        m = Model(with_optimizer(Cbc.Optimizer))
        # variables
        @variable(m,dPd[1:nd])
        @variable(m,ndPg[1:ng])
        @variable(m,pdPg[1:ng])
        @variable(m,ug[1:ng],Bin)
        @variable(m,dTheta[1:n])
        # variable bounds
        @constraint(m,-Pd .<= dPd .<= 0)
        @constraint(m, ug.*(Pg_min-Pg) .<= ndPg)
        @constraint(m, ug.*(Pg_min-Pg) .<= ndPg)
        @constraint(m, pdPg .>= 0)
        @constraint(m, ndPg .<= 0)
        @constraint(m,dTheta[1] == 0)
        # objective
        @objective(m,Max,sum(dPd)+0.01*(sum(ndPg)-sum(pdPg)+sum(ug))) # serve as much load as possible
        # mapping matrix to map loads/gens to buses
        # Power balance equality constraint
        @constraint(m,B*dTheta .== G_bus*ndPg-G_bus*pdPg-D_bus*dPd);#M_G*dPg - M_D*dPd)
        # Power flow constraints
        @constraint(m,-flow_max .<= flow0 + Xinv.*(dTheta[F] - dTheta[T]) .<= flow_max)
        ### solve the model ###
        optimize!(m);
        # collect/return the outputs
        sol_dPd=value.(dPd)
        sol_ndPg=value.(ndPg)
        sol_pdPg=value.(pdPg)
        sol_dPg = sol_pdPg+sol_ndPg;
        dPd_star = sol_dPd.*ps.baseMVA
        dPg_star = sol_dPg.*ps.baseMVA
        ps.shunt.P += dPd_star; #changes ps structure
        ps.gen.Pg[gst] += dPg_star; #changes ps structure
    else
        if (!isempty(ps.gen[gst,:]) && isempty(ps.shunt)) || (isempty(ps.gen[gst,:]) && !isempty(ps.shunt)) || (isempty(ps.gen[gst,:]) && isempty(ps.shunt))
            ps.gen.Pg[gst]  .= ps.gen.Pg[gst,:].*0.0;
            ps.shunt.P .= ps.shunt.P.*0.0;
        elseif !isempty(ps.gen) && !isempty(ps.shunt)
            Pd = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
            Pg = ps.gen.Pg[gst] ./ ps.baseMVA .* ps.gen.status[gst]
            Pg_cap = ps.gen.Pmax[gst] ./ ps.baseMVA .* ps.gen.status[gst]
            if sum(Pg_cap) >= sum(Pd)
                deltaPd = 0.0;
                deltaPg = sum(Pd)-sum(Pg);
            else
                deltaPd = sum(Pg_cap)-sum(Pd);
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
            ps.shunt.P .+= deltaPd_star;
            ps.gen.Pg[gst]  .+= deltaPg_star;
        else
            ps.shunt.P = ps.shunt.P.*0.0;
            ps.gen.Pg[gst]  = ps.gen.Pg.*0.0;
        end
    end
    #adding criteria that should produce errors if incorrect.
    @assert abs(sum(ps.shunt.P)-sum(ps.gen.Pg[gst]))<=2*tolerance
    #@assert 0.<=ps.shunt.P
    return ps
end
#end