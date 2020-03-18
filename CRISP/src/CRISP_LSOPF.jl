using JuMP
using Clp
using Cbc
using Gurobi
using SparseArrays
using LinearAlgebra

#export run_dcpf
function crisp_dcpf!(ps)
    # constants
    tolerance = 1e-4
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
    ng = size(ps.gen,1)
    G = bi[ps.gen.bus]
    Pg = ps.gen.Pg ./ ps.baseMVA .* ps.gen.status
    if any(G.<1) || any(G.>n)
        error("Bad indices in gen matrix")
    end
    Pg_bus = Array(sparse(G,ones(size(G)),Pg,n,1))
    # branch data
    brst = (ps.branch.status.==1)
    F = bi[ps.branch[brst,:f]]
    T = bi[ps.branch[brst,:t]]
    Xinv = (1 ./ ps.branch.X[brst])
    Bdc = sparse(F,T,-Xinv,n,n) + sparse(T,F,-Xinv,n,n) +
          sparse(T,T,+Xinv,n,n) + sparse(F,F,+Xinv,n,n)
    #find reference bus
    if isempty(ps.bus.id[ps.bus.bus_type.==3]) # note in Pavan's ps structure I beleive it's called 'kind' not 'bus_type'
        if isempty(ps.gen) || isempty(ps.shunt)
            ps.branch.Pf[brst] .= 0
            ps.branch.Pt[brst] .= 0
            ps.branch.Qf[brst] .= 0
            ps.branch.Qt[brst] .= 0
            if isempty(ps.shunt) && !isempty(ps.gen)
                ps.gen.Pg .= 0.0
            elseif !isempty(ps.shunt) && isempty(ps.gen)
                ps.shunt.P .= 0.0
            end
            return ps #TODO
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
        is_refgen = (ps.gen.bus.==refbusid)
        if sum(is_refgen) != 1
            println("Multiple Gen on ref bus, splitting mismatch among them.")
        #if sum(is_refgen) != 1
        #    error("Must be exactly one ref generator")
        #end
            ps.gen.Pg[is_refgen] .-= (mismatch.*ps.baseMVA)/sum(is_refgen)
        else
            ps.gen.Pg[is_refgen] .-= (mismatch.*ps.baseMVA)
        end
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

##CRISP_LSOPF1.jl
function crisp_lsopf!(ps)
    # constants
    tolerance = 1e-4
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    if n>1
        bi = sparse(ps.bus.id,fill(1,n),collect(1:n)) # helps us to find things
        # load data
        nd = size(ps.shunt,1)
        D = bi[ps.shunt.bus]
        D_bus = sparse(D,collect(1:nd),1.,n,nd);
        Pd = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
        if any(D.<1) || any(D.>n)
            error("Bad indices in shunt matrix")
        end
        #Pd_bus = Array(sparse(D,ones(size(D)),Pd,n,1))
        # gen data
        ng = size(ps.gen,1)
        G = bi[ps.gen.bus]
        G_bus = sparse(G,collect(1:ng),1.,n,ng);
        Pg = ps.gen.Pg ./ ps.baseMVA .* ps.gen.status
        Pg_max = ps.gen.Pmax ./ ps.baseMVA .* ps.gen.status
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
        m = Model(with_optimizer(Gurobi.Optimizer))
        # variables
        @variable(m,dPd[1:nd])
        @variable(m,dPg[1:ng])
        @variable(m,dTheta[1:n])
        # variable bounds
        @constraint(m,-Pd.<=dPd.<=0)
        @constraint(m,-Pg.<=dPg.<=Pg_max-Pg)
        @constraint(m,dTheta[1] == 0)
        # objective
        @objective(m,Max,sum(dPd)) # serve as much load as possible
        # mapping matrix to map loads/gens to buses
        # Power balance equality constraint
        @constraint(m,B*dTheta .== G_bus*dPg-D_bus*dPd);#M_G*dPg - M_D*dPd)
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
        bi = sparse(ps.bus.id,fill(1,n),collect(1:n));
        # load data
        nd = size(ps.shunt,1)
        D = bi[ps.shunt.bus]
        D_bus = sparse(D,collect(1:nd),1.,n,nd);
        Pd1 = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
        if any(D.<1) || any(D.>n)
            error("Bad indices in shunt matrix")
        end
        # gen data
        ng = size(ps.gen.Pg,1)
        G = bi[ps.gen.bus]
        G_bus = sparse(G,collect(1:ng),1.,n,ng);
        Pg1 = ps.gen.Pg ./ ps.baseMVA .* ps.gen.status
        ug1 = ones(ng); ug1[Pg1.==0] .= 0;
        #Pg = ps.gen[gst,:Pg] ./ ps.baseMVA .* ps.gen[gst,:status]
        Pg_max = ps.gen.Pmax ./ ps.baseMVA .* ps.gen.status
        Pg_min = ps.gen.Pmin ./ ps.baseMVA .* ps.gen.status
        if any(G.<1) || any(G.>n)
            error("Bad indices in gen matrix")
        end
        m = Model(with_optimizer(Gurobi.Optimizer))
        # variables
        @variable(m, Pd[1:nd]) # demand
        @variable(m, Pg[1:ng]) # generation
        # variable bounds constraints
        @constraint(m, 0.0 .<= Pd .<= ps.shunt.P./ps.baseMVA) # load served limits
        @constraint(m, 0 .<= Pg .<= Pg_max) # generator power limits
        # power balance
        @constraint(m, 0 .== G_bus*Pg-D_bus*Pd)
        # objective
        @objective(m, Max, 100*sum(sum(Pd)));
        ## SOLVE! ##
        optimize!(m)
        sol_Pd=value.(Pd)
        sol_Pg=value.(Pg)
        dPd_star = sol_Pd.*ps.baseMVA
        dPg_star = sol_Pg.*ps.baseMVA
        ps.shunt.P .== dPd_star; #changes ps structure
        ps.gen.Pg .== dPg_star; #changes ps structure
    end
    #adding criteria that should produce errors if incorrect.
    @assert abs(sum(ps.shunt.P)-sum(ps.gen.Pg))<=2*tolerance
    #@assert 0.<=ps.shunt.P
    return ps
end


## SECOND ORDER CONE relaxation
#using Ipopt
#=include("../src/parser.jl")
# load the case data
ps = mp2ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case6ww.m") #("../data/case6ww.m")
#ps = mp2ps("../data/case6ww.m")
ps.branch[2,:status]=0;
ps.branch[5,:status]=0; =#
function crisp_soc_lsopf(ps)
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    bi = sparse(ps.bus.id,fill(1,n),collect(1:n)) # helps us to find things
    vup = ps.bus.Vmax;#upper voltage bounds
    vlo = ps.bus.Vmin;#lower voltage bounds
    # load data
    nd = size(ps.shunt,1)
    D = bi[ps.shunt.bus]
    Pd = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
    Qd = ps.shunt.Q ./ ps.baseMVA .* ps.shunt.status
    # gen data
    ng = size(ps.gen,1)
    G = bi[ps.gen.bus]
    Pg = ps.gen.Pg ./ ps.baseMVA .* ps.gen.status
    Pmax = ps.gen.Pmax ./ ps.baseMVA .* ps.gen.status
    Pmin = ps.gen.Pmin ./ ps.baseMVA .* ps.gen.status
    Qmax = ps.gen.Qmax ./ ps.baseMVA .* ps.gen.status
    Qmin = ps.gen.Qmin ./ ps.baseMVA .* ps.gen.status
    # branch data
    brst = (ps.branch.status.==1)
    F = bi[ps.branch.f[brst]]
    T = bi[ps.branch.t[brst]]
    nl = length(T);
    flow0 = ps.branch.Pf[brst]./ps.baseMVA
    flow_max = ps.branch.rateA[brst]./ps.baseMVA # this could also be rateB
    Xinv = (1 ./ ps.branch.X[brst])
    Zinv = (1 ./ (ps.branch.R[brst] + im.*ps.branch.X[brst]))
    B = sparse(F,T,-Xinv,n,n) +
        sparse(T,F,-Xinv,n,n) +
        sparse(T,T,+Xinv,n,n) +
        sparse(F,F,+Xinv,n,n)
    bc = +ps.branch.B[brst];
    # Admittance matrix
    Yij = sparse(F,T,-Zinv, n,n) +
          sparse(T,F,-Zinv, n,n) +
          sparse(T,T,+Zinv,n,n) +
          sparse(F,F,+Zinv,n,n)
    bij = real(-Zinv);
    gij = imag(-Zinv);
    Flows = sparse(F,T,1, n,n) +
        sparse(T,F,1, n,n);
    Ys = ps.bus.Gs + im.*ps.bus.Bs# line charging admittance matrix.
    ### Build the optimization model ###
    m = Model(with_optimizer(Gurobi.Optimizer))
    #m = Model(solver = IpoptSolver())
    # variables
    @variable(m,Pg[1:ng])
    @variable(m,Qg[1:ng])
    @variable(m,dTheta[1:n])
    @variable(m,Pij[1:nl])
    @variable(m,Qij[1:nl])
    @variable(m,Pji[1:nl])
    @variable(m,Qji[1:nl])
    @variable(m,W_ii[1:n])
    @variable(m,Wr_ij[1:nl])
    @variable(m,Wi_ij[1:nl])
    @variable(m,Ws_ii[1:n])
    @variable(m, z_v[1:n])
    @variable(m, z_g[1:ng])
    @variable(m, z_d[1:nd])
    @variable(m, z_s[1:n])
    # variable bounds
    @constraint(m, Pg .<= Pmax)
    @constraint(m, Pg .>= 0)
    @constraint(m, Qg .<= Qmax)
    @constraint(m, Qg .>= Qmin)
    @constraint(m, 0 .<= z_v)
    @constraint(m, 0 .<= z_g)
    @constraint(m, 0 .<= z_d)
    @constraint(m, 0 .<= z_s)
    @constraint(m, 1 .>= z_v)
    @constraint(m, 1 .>= z_g)
    @constraint(m, 1 .>= z_d)
    @constraint(m, 1 .>= z_s)
#    @constraint(m,dTheta[1] == 0)
    @constraint(m, constr1[i=1:n], W_ii[i] <= z_v[i]*vup[i]^2)
    @constraint(m, constr2[i=1:n], W_ii[i] >= z_v[i]*vlo[i]^2)
    #@constraint(m, (Wr_ij.^2+Wi_ij.^2) .<= W_ii[F].*W_ii[T])
    # McCormick envelope representing Ws_ii = z_s.*Wii
    #@constraint(m, constr[i=1:n], Ws_ii[i] >= z_s[i]*vlo[i]^2)
    #@constraint(m, constr[i=1:n], Ws_ii[i] >= W_ii[i] + z_s[i]*vup[i]^2 - vup[i]^2)
    #@constraint(m, constr[i=1:n], Ws_ii[i] <= z_s[i]*vup[i]^2)
    #@constraint(m, constr[i=1:n], Ws_ii[i] <= W_ii[i] + z_s[i]*vlo[i]^2 - vlo[i]^2)
    # objective: serve as much load as possible and keep as much of the network connected as possible
    C = 100;
    @objective(m,Min,-sum(z_v)-sum(z_g)-sum(z_s)-C*sum(z_d.*Pd))
    # mapping matrix to map loads/gens to buses
    M_D = sparse(D,1:nd,1.0,n,nd)
    M_G = sparse(G,1:ng,1.0,n,ng)
    # Power balance equality constraint
    Ysr = real(Ys);
    Ysi = imag(Ys);
    @constraint(m, constr3[i=1:n], (M_G*Pg)[i] - (M_D*z_d)[i].*(M_D*Pd)[i]
                     - Ysr[i].*W_ii[i] ==  sum(Pij[F.==i]) - sum(Pji[T.==i]))
    @constraint(m, constr4[i=1:n], (M_G*Qg)[i] - (M_D*z_d)[i].*(M_D*Qd)[i]
                    - Ysi[i].*W_ii[i] ==  sum(Qij[F.==i]) - sum(Qji[T.==i])) #
    # Power flow constraints
    @constraint(m, Pij.^2+Qij.^2 .<= flow_max.^2)
    @constraint(m, Pji.^2+Qji.^2 .<= flow_max.^2)
    @constraint(m, Pij .== (bij).*W_ii[F] - (bij.*Wr_ij+gij.*Wi_ij))
    @constraint(m, Pji .== (bij).*W_ii[T] + (bij.*Wr_ij-gij.*Wi_ij))
    @constraint(m, Qij .== (-gij - bc/2).*W_ii[F] + (gij.*Wr_ij-bij.*Wi_ij))
    @constraint(m, Qji .== (-gij - bc/2).*W_ii[T] + (gij.*Wr_ij+bij.*Wi_ij))
    @constraint(m,tan(-pi/2).*Wr_ij .<= Wi_ij)
    @constraint(m,tan(pi/2).*Wr_ij .>= Wi_ij)
    ### solve the model ###
    optimize!(m)
    # collect/return the outputs
    dPd_star = value.(z_d).*Pd.*ps.baseMVA
    dPg_star = value.(z_g).*Pg.*ps.baseMVA
    return (dPd_star, dPg_star)
end

## Added generators
## CRISP_LSOPF_gen.jl
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
    brst = (ps.branch.status.==1)
    F = bi[ps.branch.f[brst]]
    T = bi[ps.branch.t[brst]]
    Xinv = (1 ./ ps.branch.X[brst])
    Bdc = sparse(F,T,-Xinv,n,n) + sparse(T,F,-Xinv,n,n) +
          sparse(T,T,+Xinv,n,n) + sparse(F,F,+Xinv,n,n)
    #find reference bus
    if isempty(ps.bus.id[ps.bus.bus_type.==3]) # note in Pavan's ps structure I beleive it's called 'kind' not 'bus_type'
        if (isempty(ps.gen[gst,:]) && isempty(ps.stoarage)) || (isempty(ps.shunt) && isempty(ps.stoarage)) || (isempty(ps.shunt) && isempty(ps.gen[gst,:]))
            ps.branch.Pf[brst] .= 0
            ps.branch.Pt[brst] .= 0
            ps.branch.Qf[brst] .= 0
            ps.branch.Qt[brst] .= 0
            if !isempty(ps.gen[gst,:])
                ps.gen.Pg[gst] .= 0.0
            elseif !isempty(ps.shunt)
                ps.shunt.P .= 0.0
            elseif !isempty(ps.storage)
                ps.storage.Ps .= 0.0
            end
            return ps #TODO
        elseif !isempty(ps.gen)
            maxGen = findmax(ps.gen.Pg)[2]
            busID = ps.gen.bus[maxGen];
            isref = (busID.==ps.bus.id)
            nonref = .~isref
        end
    elseif sum(ps.gen.status[ps.gen.bus.==ps.bus.id[ps.bus.bus_type.==3]])==0
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
            maxGen = findmax(ps.gen.Pg)[2]
            busID = ps.gen[maxGen,:bus];
            isref = (busID.==ps.bus.id)
            nonref = .~isref
        end
    else
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
            isref = (ps.bus.bus_type.==3)
            nonref = .~isref
        end
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
        refbusid = ps.bus.id[isref]
        is_refgen = (ps.gen.bus.==refbusid)
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
        D = bi[ps.shunt.bus]
        D_bus = sparse(D,collect(1:nd),1.,n,nd);
        Pd = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
        if any(D.<1) || any(D.>n)
            error("Bad indices in shunt matrix")
        end
        #Pd_bus = Array(sparse(D,ones(size(D)),Pd,n,1))
        # gen data
        gst = (ps.gen.status.==1)
        ng = size(ps.gen[gst,:],1)
        G = bi[ps.gen.bus[gst]]
        G_bus = sparse(G,collect(1:ng),1.,n,ng);
        Pg = ps.gen.Pg[gst] ./ ps.baseMVA .* ps.gen.status[gst]
        Pg_max = ps.gen.Pmax[gst] ./ ps.baseMVA .* ps.gen.status[gst]
        Pg_min = ps.gen.Pmin[gst] ./ ps.baseMVA .* ps.gen.status[gst]
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
        m = Model(with_optimizer(Gurobi.Optimizer))
        # variables
        @variable(m,dPd[1:nd])
        @variable(m,ndPg[1:ng])
        @variable(m,pdPg[1:ng])
        @variable(m,ug[1:ng],Bin)
        @variable(m,dTheta[1:n])
        # variable bounds
        @constraint(m,-Pd .<= dPd .<= 0)
        @constraint(m, ug.*(Pg_min) .<= Pg+ndPg+pdPg)
        @constraint(m, ug.*(Pg_max) .>= Pg+pdPg+ndPg)
        @constraint(m, pdPg .>= 0)
        @constraint(m, ndPg .<= 0)
        @constraint(m,dTheta[1] == 0)
        # objective
        @objective(m,Max,sum(dPd)+0.01*(0.1*sum(ndPg)-0.1*sum(pdPg)+sum(ug))) # serve as much load as possible
        # mapping matrix to map loads/gens to buses
        # Power balance equality constraint
        @constraint(m,B*dTheta .== G_bus*ndPg+G_bus*pdPg-D_bus*dPd);#M_G*dPg - M_D*dPd)
        # Power flow constraints
        @constraint(m,-flow_max .<= flow0 + Xinv.*(dTheta[F] - dTheta[T]) .<= flow_max)
        ### solve the model ###
        optimize!(m);
        # collect/return the outputs
        sol_dPd=value.(dPd)
        sol_ndPg=value.(ndPg)
        sol_pdPg=value.(pdPg)
        dPd_star = sol_dPd.*ps.baseMVA
        dPg_star = (sol_pdPg+sol_ndPg).*ps.baseMVA
        ps.shunt.P += dPd_star; #changes ps structure
        ps.gen.Pg[gst] += dPg_star; #changes ps structure
    else
        gst = (ps.gen.status.==1)
        Pd = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
        nd = length(Pd)
        Pg = ps.gen.Pg[gst] ./ ps.baseMVA .* ps.gen.status[gst]
        ng = length(Pg)
        Pg_min = ps.gen.Pmin[gst] ./ ps.baseMVA .* ps.gen.status[gst]
        Pg_max = ps.gen.Pmax[gst] ./ ps.baseMVA .* ps.gen.status[gst]
        if !isempty(Pd) && !isempty(Pg)
            m = Model(with_optimizer(Gurobi.Optimizer))
            # variables
            @variable(m,dPd[1:nd])
            @variable(m,ndPg[1:ng])
            @variable(m,pdPg[1:ng])
            @variable(m,ug[1:ng],Bin)
            # variable bounds
            @constraint(m,-Pd .<= dPd .<= 0)
            @constraint(m, ug.*(Pg_min) .<= Pg+ndPg+pdPg)
            @constraint(m, ug.*(Pg_max) .>= Pg+pdPg+ndPg)
            @constraint(m, pdPg .>= 0)
            @constraint(m, ndPg .<= 0)
            # objective
            @objective(m,Max,sum(dPd)+0.01*(0.1*sum(ndPg)-0.1*sum(pdPg)+sum(ug))) # serve as much load as possible
            # mapping matrix to map loads/gens to buses
            # Power balance equality constraint
            @constraint(m,sum(Pg+ndPg+pdPg)-sum(Pd+dPd) == 0);
            ### solve the model ###
            optimize!(m);
            sol_dPd=value.(dPd)
            sol_ndPg=value.(ndPg)
            sol_pdPg=value.(pdPg)
            dPd_star = sol_dPd.*ps.baseMVA
            dPg_star = (sol_pdPg+sol_ndPg).*ps.baseMVA
            ps.shunt.P .+= dPd_star;
            ps.gen.Pg[gst]  .+= dPg_star;
        else
            ps.shunt.P = ps.shunt.P.*0.0;
            ps.gen.Pg  = ps.gen.Pg.*0.0;
        end
    end
    #adding criteria that should produce errors if incorrect.
    @assert abs(sum(ps.shunt.P)-sum(ps.gen.Pg[gst]))<=2*tolerance
    #@assert 0.<=ps.shunt.P
    return ps
end

function crisp_dcpf_g_s1!(ps)
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
    brst = (ps.branch.status.==1)
    F = bi[ps.branch.f[brst]]
    T = bi[ps.branch.t[brst]]
    Xinv = (1 ./ ps.branch.X[brst])
    Bdc = sparse(F,T,-Xinv,n,n) + sparse(T,F,-Xinv,n,n) +
          sparse(T,T,+Xinv,n,n) + sparse(F,F,+Xinv,n,n)
    #find reference bus
    if isempty(ps.bus.id[ps.bus.bus_type.==3]) # note in Pavan's ps structure I beleive it's called 'kind' not 'bus_type'
        if (isempty(ps.gen[gst,:]))  || (isempty(ps.shunt) && isempty(ps.storage))
            ps.branch.Pf[brst] .= 0
            ps.branch.Pt[brst] .= 0
            ps.branch.Qf[brst] .= 0
            ps.branch.Qt[brst] .= 0
            if !isempty(ps.gen[gst,:])
                ps.gen.Pg[gst] .= 0.0
            end
            if !isempty(ps.shunt)
                ps.shunt.P .= 0.0
            end
            if !isempty(ps.storage)
                ps.storage.Ps .= 0.0
            end
            return ps #TODO
        elseif !isempty(ps.gen[gst,:])
            maxGen = findmax(ps.gen.Pg)[2]
            busID = ps.gen[maxGen,:bus];
            isref = (busID.==ps.bus.id)
            nonref = .~isref
        end
    elseif sum(ps.gen.status[ps.gen.bus.==ps.bus.id[ps.bus.bus_type.==3]])==0
        if isempty(ps.gen[gst,:]) || (isempty(ps.shunt) && isempty(ps.storage))
            ps.branch.Pf[brst] .= 0
            ps.branch.Pt[brst] .= 0
            ps.branch.Qf[brst] .= 0
            ps.branch.Qt[brst] .= 0
            if !isempty(ps.gen[gst,:])
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
        if isempty(ps.gen[gst,:]) || (isempty(ps.shunt) && isempty(ps.storage))
            ps.branch.Pf[brst] .= 0
            ps.branch.Pt[brst] .= 0
            ps.branch.Qf[brst] .= 0
            ps.branch.Qt[brst] .= 0
            if !isempty(ps.gen[gst,:])
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
        refbusid = ps.bus[isref,:id]
        is_refgen = (ps.gen.bus.==refbusid)
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

function crisp_lsopf_g_s1!(ps,dt)
    # constants
    tolerance = 1e-6
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    if n>1
        bi = sparse(ps.bus.id,fill(1,n),collect(1:n)) # helps us to find things
        # load data
        nd = size(ps.shunt,1)
        D = bi[ps.shunt.bus]
        D_bus = sparse(D,collect(1:nd),1.,n,nd);
        Pd = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
        if any(D.<1) || any(D.>n)
            error("Bad indices in shunt matrix")
        end
        #Pd_bus = Array(sparse(D,ones(size(D)),Pd,n,1))
        # gen data
        gst = (ps.gen.status.==1)
        ng = size(ps.gen[gst,:],1)
        G = bi[ps.gen.bus[gst]]
        G_bus = sparse(G,collect(1:ng),1.,n,ng);
        Pg = ps.gen.Pg[gst] ./ ps.baseMVA .* ps.gen.status[gst]
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
            error("Bad indices in sstorage matrix")
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
        # variables
        @variable(m,dPd[1:nd])
        @variable(m,dPs[1:ns])
        @variable(m,ndPg[1:ng])
        @variable(m,pdPg[1:ng])
        @variable(m,ug[1:ng],Bin)
        @variable(m,dTheta[1:n])
        # variable bounds
        @constraint(m,-Pd .<= dPd .<= 0)
        @constraint(m, Ps_min .<= Ps+dPs .<= Ps_max)
        @constraint(m, 0.01 .<= E + (Ps+dPs).*dt .<= E_max)
        @constraint(m, ug.*(Pg_min) .<= Pg+ndPg+pdPg)
        @constraint(m, ug.*(Pg_max) .>= Pg+pdPg+ndPg)
        @constraint(m, pdPg .>= 0)
        @constraint(m, ndPg .<= 0)
        @constraint(m,dTheta[1] == 0)
        # objective
        @objective(m,Max,sum(dPd)+0.01*(0.1*sum(ndPg)-0.1*sum(pdPg)+sum(ug))) # serve as much load as possible
        # mapping matrix to map loads/gens to buses
        # Power balance equality constraint
        @constraint(m,B*dTheta .== G_bus*ndPg+G_bus*pdPg+S_bus*dPs-D_bus*dPd);#M_G*dPg - M_D*dPd)
        # Power flow constraints
        @constraint(m,-flow_max .<= flow0 + Xinv.*(dTheta[F] - dTheta[T]) .<= flow_max)
        ### solve the model ###
        optimize!(m);
        # collect/return the outputs
        sol_dPd=value.(dPd)
        sol_dPs=value.(dPs)
        sol_ndPg=value.(ndPg)
        sol_pdPg=value.(pdPg)
        dPd_star = sol_dPd.*ps.baseMVA
        dPs_star = sol_dPs.*ps.baseMVA
        dPg_star = (sol_pdPg+sol_ndPg).*ps.baseMVA
        ps.shunt.P += dPd_star; #changes ps structure
        ps.storage.Ps += dPs_star;
        ps.storage.E += ps.storage.Ps.*dt;
        ps.gen.Pg[gst] += dPg_star; #changes ps structure
    else
        gst = (ps.gen.status.==1)
        Pd = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
        nd = length(Pd)
        Pg = ps.gen.Pg[gst] ./ ps.baseMVA .* ps.gen.status[gst]
        ng = length(Pg)
        Pg_min = ps.gen.Pmin[gst] ./ ps.baseMVA .* ps.gen.status[gst]
        Pg_max = ps.gen.Pmax[gst] ./ ps.baseMVA .* ps.gen.status[gst]
        Ps = ps.storage.Ps ./ ps.baseMVA .* ps.storage.status
        ns = length(Ps);
        E = ps.storage.E ./ ps.baseMVA .* ps.storage.status
        E_max = ps.storage.Emax ./ ps.baseMVA .* ps.storage.status
        Ps_max = ps.storage.Ps ./ ps.baseMVA .* ps.storage.status
        Ps_min = ps.storage.Ps ./ ps.baseMVA .* ps.storage.status
        if sum(isempty(Pd) + isempty(Pg) + isempty(Ps)) <=1
            m = Model(with_optimizer(Gurobi.Optimizer))
            # variables
            @variable(m,dPd[1:nd])
            @variable(m,dPs[1:ns])
            @variable(m,ndPg[1:ng])
            @variable(m,pdPg[1:ng])
            @variable(m,ug[1:ng],Bin)
            # variable bounds
            @constraint(m,-Pd .<= dPd .<= 0)
            @constraint(m, Ps_min .<= Ps+dPs .<= Ps_max)
            @constraint(m, 0 .<= E + (Ps+dPs).*dt .<= E_max)
            @constraint(m, ug.*(Pg_min) .<= Pg+ndPg+pdPg)
            @constraint(m, ug.*(Pg_max) .>= Pg+pdPg+ndPg)
            @constraint(m, pdPg .>= 0)
            @constraint(m, ndPg .<= 0)
            # objective
            @objective(m,Max,sum(dPd)+0.01*(0.1*sum(ndPg)-0.1*sum(pdPg)+sum(ug))) # serve as much load as possible
            # mapping matrix to map loads/gens to buses
            # Power balance equality constraint
            @constraint(m,sum(Pg+ndPg+pdPg)+sum(Ps+dPs)-sum(Pd+dPd) == 0);
            ### solve the model ###
            optimize!(m);
            sol_dPd=value.(dPd)
            sol_dPs=value.(dPs)
            sol_ndPg=value.(ndPg)
            sol_pdPg=value.(pdPg)
            dPd_star = sol_dPd.*ps.baseMVA
            dPs_star = sol_dPs.*ps.baseMVA
            dPg_star = (sol_pdPg+sol_ndPg).*ps.baseMVA
            ps.shunt.P .+= dPd_star;
            ps.storage.Ps += dPs_star;
            ps.storage.E += ps.storage.Ps.*dt;
            ps.gen.Pg[gst]  .+= dPg_star;
        else
            ps.storage.Ps = ps.storage.Ps.*0.0;
            ps.shunt.P = ps.shunt.P.*0.0;
            ps.gen.Pg  = ps.gen.Pg.*0.0;
        end
    end
    #adding criteria that should produce errors if incorrect.
    @assert abs(sum(ps.shunt.P)-sum(ps.storage.Ps)-sum(ps.gen.Pg[gst]))<=2*tolerance
    #@assert sum(ps.storage.E .< 0)==0
    #@assert 0.<=ps.shunt.P
    return ps
end
#end



## ADDED GENERATORS AND storage
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
        if (isempty(ps.gen[gst]))  || (isempty(ps.shunt) && isempty(ps.storage))
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
        elseif !isempty(ps.gen[gst])
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
        if isempty(ps.gen[gst,:]) || (isempty(ps.shunt) && isempty(ps.storage))
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

function crisp_lsopf_g_s1!(ps,dt)
    # constants
    tolerance = 1e-6
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    if n>1
        bi = sparse(ps.bus.id,fill(1,n),collect(1:n)) # helps us to find things
        # load data
        nd = size(ps.shunt,1)
        D = bi[ps.shunt.bus]
        D_bus = sparse(D,collect(1:nd),1.,n,nd);
        Pd = ps.shunt.P  ./ ps.baseMVA .* ps.shunt.status
        if any(D.<1) || any(D.>n)
            error("Bad indices in shunt matrix")
        end
        #Pd_bus = Array(sparse(D,ones(size(D)),Pd,n,1))
        # gen data
        gst = (ps.gen.status.==1)
        ng = size(ps.gen[gst],1)
        G = bi[ps.gen.bus[gst]]
        G_bus = sparse(G,collect(1:ng),1.,n,ng);
        Pg = ps.gen.Pg[gst] ./ ps.baseMVA .* ps.gen.status[gst]
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
            error("Bad indices in sstorage matrix")
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
        # variables
        @variable(m,dPd[1:nd])
        @variable(m,dPs[1:ns])
        @variable(m,ndPg[1:ng])
        @variable(m,pdPg[1:ng])
        @variable(m,ug[1:ng],Bin)
        @variable(m,dTheta[1:n])
        # variable bounds
        @constraint(m,-Pd .<= dPd .<= 0)
        @constraint(m, Ps_min .<= Ps+dPs .<= Ps_max)
        @constraint(m, 0.01 .<= E + (Ps+dPs).*dt/60 .<= E_max)
        @constraint(m, ug.*(Pg_min) .<= Pg+ndPg+pdPg)
        @constraint(m, ug.*(Pg_max) .>= Pg+pdPg+ndPg)
        @constraint(m, pdPg .>= 0)
        @constraint(m, ndPg .<= 0)
        @constraint(m,dTheta[1] == 0)
        # objective
        @objective(m,Max,sum(dPd)+0.01*(0.1*sum(ndPg)-0.1*sum(pdPg)+sum(ug))) # serve as much load as possible
        # mapping matrix to map loads/gens to buses
        # Power balance equality constraint
        @constraint(m,B*dTheta .== G_bus*ndPg+G_bus*pdPg+S_bus*dPs-D_bus*dPd);#M_G*dPg - M_D*dPd)
        # Power flow constraints
        @constraint(m,-flow_max .<= flow0 + Xinv.*(dTheta[F] - dTheta[T]) .<= flow_max)
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
    else
        gst = (ps.gen.status .== 1)
        Pd = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
        nd = length(Pd)
        Pg = ps.gen.Pg[gst] ./ ps.baseMVA .* ps.gen.status[gst]
        ng = length(Pg)
        Pg_min = ps.gen.Pmin[gst] ./ ps.baseMVA .* ps.gen.status[gst]
        Pg_max = ps.gen.Pmax[gst] ./ ps.baseMVA .* ps.gen.status[gst]
        Ps = ps.storage.Ps  ./ ps.baseMVA .* ps.storage.status
        ns = length(Ps);
        E = ps.storage.E  ./ ps.baseMVA .* ps.storage.status
        E_max = ps.storage.Emax  ./ ps.baseMVA .* ps.storage.status
        Ps_max = ps.storage.Ps  ./ ps.baseMVA .* ps.storage.status
        Ps_min = ps.storage.Ps  ./ ps.baseMVA .* ps.storage.status
        if sum(isempty(Pd) + isempty(Pg) + isempty(Ps)) <=1
            m = Model(with_optimizer(Gurobi.Optimizer))
            # variables
            @variable(m,dPd[1:nd])
            @variable(m,dPs[1:ns])
            @variable(m,ndPg[1:ng])
            @variable(m,pdPg[1:ng])
            @variable(m,ug[1:ng],Bin)
            # variable bounds
            @constraint(m,-Pd .<= dPd .<= 0)
            @constraint(m, Ps_min .<= Ps+dPs .<= Ps_max)
            @constraint(m, 0 .<= E + (Ps+dPs).*dt .<= E_max)
            @constraint(m, ug.*(Pg_min) .<= Pg+ndPg+pdPg)
            @constraint(m, ug.*(Pg_max) .>= Pg+pdPg+ndPg)
            @constraint(m, pdPg .>= 0)
            @constraint(m, ndPg .<= 0)
            # objective
            @objective(m,Max,sum(dPd)+0.01*(0.1*sum(ndPg)-0.1*sum(pdPg)+sum(ug))) # serve as much load as possible
            # mapping matrix to map loads/gens to buses
            # Power balance equality constraint
            @constraint(m,sum(Pg+ndPg+pdPg)+sum(Ps+dPs)-sum(Pd+dPd) == 0);
            ### solve the model ###
            optimize!(m);
            sol_dPd=value.(dPd)
            sol_dPs=value.(dPs)
            sol_ndPg=value.(ndPg)
            sol_pdPg=value.(pdPg)
            dPd_star = sol_dPd.*ps.baseMVA
            dPs_star = sol_dPs.*ps.baseMVA
            dPg_star = (sol_pdPg+sol_ndPg).*ps.baseMVA
            ps.shunt.P .+= dPd_star;
            ps.storage.Ps += dPs_star;
            ps.storage.E += ps.storage.Ps.*dt;
            ps.gen.Pg[gst]  .+= dPg_star;
        else
            ps.storage.Ps = ps.storage.Ps.*0.0;
            ps.shunt.P = ps.shunt.P.*0.0;
            ps.gen.Pg  = ps.gen.Pg.*0.0;
        end
    end
    #adding criteria that should produce errors if incorrect.
    @assert abs(sum(ps.shunt.P)-sum(ps.storage.Ps)-sum(ps.gen.Pg[gst]))<=2*tolerance
    #@assert sum(ps.storage.E .< 0)==0
    #@assert 0.<=ps.shunt.P
    return ps
end
#end
function crisp_lsopf_g_s!(ps,dt)
t = [0; dt];
Ti = size(t,1);
# constants
tolerance = 1e-4
### collect the data that we will need ###
# bus data
n = size(ps.bus,1) # the number of buses
if n>1
    bi = ps.bi;
    # load data
    nd = size(ps.shunt,1)
    D = bi[ps.shunt.bus]
    D_bus = sparse(D,collect(1:nd),1.,n,nd);
    Pd1 = ps.shunt.P  ./ ps.baseMVA .* ps.shunt.status
    if any(D.<1) || any(D.>n)
        error("Bad indices in shunt matrix")
    end
    # gen data
    gst = (ps.gen.status .== 1);
    ng = size(ps.gen.Pg[gst],1)
    G = bi[ps.gen.bus[gst]]
    G_bus = sparse(G,collect(1:ng),1.,n,ng);
    Pg1 = ps.gen.Pg[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    ug1 = ones(ng); ug1[Pg1.==0] .= 0;
    #Pg = ps.gen[gst,:Pg] ./ ps.baseMVA .* ps.gen[gst,:status]
    Pg_max = ps.gen.Pmax[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    Pg_min = ps.gen.Pmin[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    if any(G.<1) || any(G.>n)
        error("Bad indices in gen matrix")
    end
    # storage data
    ns = size(ps.storage,1)
    S = bi[ps.storage.bus];
    S_bus = sparse(S,collect(1:ns),1.,n,ns);
    Ps1 = ps.storage.Ps ./ ps.baseMVA .* ps.storage.status
    E1 = ps.storage.E ./ ps.baseMVA .* ps.storage.status
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
    flow_max = ps.branch.rateA[brst]./ps.baseMVA # this could also be rateB
    Xinv = (1 ./ ps.branch.X[brst])
    B = sparse(F,T,-Xinv,n,n) +
        sparse(T,F,-Xinv,n,n) +
        sparse(T,T,+Xinv,n,n) +
        sparse(F,F,+Xinv,n,n)
    m = Model(with_optimizer(Gurobi.Optimizer))
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
    @constraint(m, E .== (E1 + ((dt/60) .*(Ps)))) # storage energy at next time step
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
    dE_star = sol_E.*ps.baseMVA
    dPd_star = sol_Pd.*ps.baseMVA./ps.shunt.P # % load served
    dPs_star = sol_Ps.*ps.baseMVA
    dPg_star = sol_Pg.*ps.baseMVA
else
    bi = ps.bi;
    # load data
    nd = size(ps.shunt,1)
    D = bi[ps.shunt.bus]
    D_bus = sparse(D,collect(1:nd),1.,n,nd);
    Pd1 = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
    if any(D.<1) || any(D.>n)
        error("Bad indices in shunt matrix")
    end
    # gen data
    gst = (ps.gen.status .== 1);
    ng = size(ps.gen.Pg[gst],1)
    G = bi[ps.gen.bus[gst]]
    G_bus = sparse(G,collect(1:ng),1.,n,ng);
    Pg1 = ps.gen.Pg[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    ug1 = ones(ng); ug1[Pg1.==0] .= 0;
    #Pg = ps.gen[gst,:Pg] ./ ps.baseMVA .* ps.gen[gst,:status]
    Pg_max = ps.gen.Pmax[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    Pg_min = ps.gen.Pmin[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    if any(G.<1) || any(G.>n)
        error("Bad indices in gen matrix")
    end
    # storage data
    ns = size(ps.storage,1)
    S = bi[ps.storage.bus];
    S_bus = sparse(S,collect(1:ns),1.,n,ns);
    Ps1 = ps.storage.Ps ./ ps.baseMVA .* ps.storage.status
    E1 = ps.storage.E ./ ps.baseMVA .* ps.storage.status
    E_max = ps.storage.Emax ./ ps.baseMVA .* ps.storage.status
    Ps_max = ps.storage.Psmax ./ ps.baseMVA .* ps.storage.status
    Ps_min = ps.storage.Psmin ./ ps.baseMVA .* ps.storage.status
    if any(S.<1) || any(S.>n)
        error("Bad indices in storage matrix")
    end
    m = Model(with_optimizer(Gurobi.Optimizer))
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

#from CRISP_LSOPF_gen_stor.jl

function crisp_lsopf_g_s2!(ps,dt) # changed function name to avoid multiple functions with the same name
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
    D = bi[ps.shunt.bus]
    D_bus = sparse(D,collect(1:nd),1.,n,nd);
    Pd1 = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
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
    ng = size(ps.gen.Pg[gst],1)
    G = bi[ps.gen.bus[gst]]
    G_bus = sparse(G,collect(1:ng),1.,n,ng);
    Pg1 = ps.gen.Pg[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    ug1 = ones(ng); ug1[Pg1.==0] .= 0;
    #Pg = ps.gen[gst,:Pg] ./ ps.baseMVA .* ps.gen[gst,:status]
    Pg_max = ps.gen.Pmax[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    Pg_min = ps.gen.Pmin[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    if any(G.<1) || any(G.>n)
        error("Bad indices in gen matrix")
    end
    # storage data
    ns = size(ps.storage,1)
    S = bi[ps.storage.bus];
    S_bus = sparse(S,collect(1:ns),1.,n,ns);
    Ps1 = ps.storage.Ps ./ ps.baseMVA .* ps.storage.status
    E1 = ps.storage.E ./ ps.baseMVA .* ps.storage.status
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
    D = bi[ps.shunt.bus]
    D_bus = sparse(D,collect(1:nd),1.,n,nd);
    Pd1 = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
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
    ng = size(ps.gen.Pg[gst],1)
    G = bi[ps.gen.bus[gst]]
    G_bus = sparse(G,collect(1:ng),1.,n,ng);
    Pg1 = ps.gen.Pg[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    Pg_max = ps.gen.Pmax[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    Pg_min = ps.gen.Pmin[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    if any(G.<1) || any(G.>n)
        error("Bad indices in gen matrix")
    end
    # storage data
    ns = size(ps.storage,1)
    S = bi[ps.storage.bus];
    S_bus = sparse(S,collect(1:ns),1.,n,ns);
    Ps1 = ps.storage.Ps ./ ps.baseMVA .* ps.storage.status
    E1 = ps.storage.E ./ ps.baseMVA .* ps.storage.status
    E_max = ps.storage.Emax ./ ps.baseMVA .* ps.storage.status
    Ps_max = ps.storage.Psmax ./ ps.baseMVA .* ps.storage.status
    Ps_min = ps.storage.Psmin ./ ps.baseMVA .* ps.storage.status
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
    D = bi[ps.shunt.bus]
    D_bus = sparse(D,collect(1:nd),1.,n,nd);
    Pd1 = ps.shunt.P ./ ps.baseMVA;
    if any(D.<1) || any(D.>n)
        error("Bad indices in shunt matrix")
    end
    # gen data
    gst = (ps.gen.status .== 1);
    ng = size(ps.gen.Pg[gst],1)
    G = bi[ps.gen.bus[gst]]
    G_bus = sparse(G,collect(1:ng),1.,n,ng);
    Pg1 = ps.gen.Pg[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    ug1 = ones(ng); ug1[Pg1.==0] .= 0;
    #Pg = ps.gen[gst,:Pg] ./ ps.baseMVA .* ps.gen[gst,:status]
    Pg_max = ps.gen.Pmax[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    Pg_min = ps.gen.Pmin[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    if any(G.<1) || any(G.>n)
        error("Bad indices in gen matrix")
    end
    # storage data
    ns = size(ps.storage,1)
    S = bi[ps.storage.bus];
    S_bus = sparse(S,collect(1:ns),1.,n,ns);
    Ps1 = ps.storage.Ps ./ ps.baseMVA .* ps.storage.status
    E1 = ps.storage.E ./ ps.baseMVA .* ps.storage.status
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
    D = bi[ps.shunt.bus]
    D_bus = sparse(D,collect(1:nd),1.,n,nd);
    Pd1 = ps.shunt.P ./ ps.baseMVA .* ps.shunt.status
    if any(D.<1) || any(D.>n)
        error("Bad indices in shunt matrix")
    end
    # gen data
    g1 = (ps.gen.status .== 1);
    g2 = (ps.gen.Pg .!=0);
    g3 = (ps.gen.time_on .== 0)
    ps.gen.time_off[!g2 && g3] += dt./60;
    gst = g1 && g2;
    ng = size(ps.gen.Pg[gst],1)
    G = bi[ps.gen.bus[gst]]
    G_bus = sparse(G,collect(1:ng),1.,n,ng);
    Pg1 = ps.gen.Pg[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    ug1 = ones(ng); ug1[Pg1.==0] .= 0;
    #Pg = ps.gen[gst,:Pg] ./ ps.baseMVA .* ps.gen[gst,:status]
    Pg_max = ps.gen.Pmax[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    Pg_min = ps.gen.Pmin[gst] ./ ps.baseMVA .* ps.gen.status[gst]
    if any(G.<1) || any(G.>n)
        error("Bad indices in gen matrix")
    end
    # storage data
    ns = size(ps.storage,1)
    S = bi[ps.storage.bus];
    S_bus = sparse(S,collect(1:ns),1.,n,ns);
    Ps1 = ps.storage.Ps ./ ps.baseMVA .* ps.storage.status
    E1 = ps.storage.E ./ ps.baseMVA .* ps.storage.status
    E_max = ps.storage.Emax./ ps.baseMVA .* ps.storage.status
    Ps_max = ps.storage.Psmax ./ ps.baseMVA .* ps.storage.status
    Ps_min = ps.storage.Psmin ./ ps.baseMVA .* ps.storage.status
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

##CRISP_LSOPF_gen1.jl
function crisp_dcpf_g1_s!(ps)
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
            busID = ps.gen.bus[maxGen];
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

function crisp_mh_lsopf_var!(ps,dt,ug,ul,Pd_max,Pg_max1,load_shed_cost;t_win=dt,w_g=0.1)
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
    Pg_max = (Pg_max1 ./ ps.baseMVA) .* ps.gen.status
    if sum(names(ps.gen) .== :RampRateMWMin) .!=0
        RR = ps.gen.RampRateMWMin .* dt ./ ps.baseMVA
    else
        RR = ps.gen.ramp30 ./30 .* dt ./ ps.baseMVA
    end
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
    m = Model(with_optimizer(Gurobi.Optimizer))
    #m = Model(with_optimizer(Cbc.Optimizer))
    # variables
    @variable(m, Pd[1:nd]) # demand
    @variable(m, Pg[1:ng]) # generation
    @variable(m, Ps[1:ns]) # power flow into or out of storage (negative flow = charging)
    @variable(m, E[1:ns, 1:2]) # energy level in battery
    # fix battery starting charge
    for s in 1:ns
        fix(E[s,1], E1[s], force = true)
    end
    # variable bounds constraints
    @constraint(m, stPdcon, 0.0 .<= Pd[:] .<= Pdmax[:,1]) # load served limits
    @constraint(m, stPscon, Ps_min .<= Ps[:] .<= Ps_max) # storage power flow
    @constraint(m, stEPscon, E[:,2] .== (E[:,1] - ((dt/60) .* (Ps[:])))) # storage energy at next time step
    @constraint(m, stEcon, min_bat_E .<= (E[:,2]) .<= E_max) # storage energy
    @constraint(m, genPgucon, Pg[:] .<= ug[:,1] .* Pg_max[:,1]) # generator power limits upper
    @constraint(m, genPglcon, 0 .<= Pg[:]) # generator power limits lower
    if n > 1
        @variable(m, Theta[1:n])
        @constraint(m, Theta[1] .== 0); # set first bus as reference bus: V angle to 0
        brst = falses(nb);
        for b in 1:nb
            if ul[b,1] == 1
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
        @constraint(m, B*Theta[:] .== G_bus*Pg[:]+S_bus*Ps[:]-D_bus*Pd[:])
        # power flow limits
        @constraint(m, -flow_max[brst] .<= Xinv .* (Theta[F] - Theta[T]) .<= flow_max[brst])
    else
        @constraint(m, PBnoTheta, 0.0 .== G_bus*Pg[:]+S_bus*Ps[:]-D_bus*Pd[:])
    end
    # objective
    @objective(m, Max, sum(Pd) + sum(w_g.*Pg)) #TODO: constants
    ## SOLVE! ##
    optimize!(m)
    sol_Pd=value.(Pd)[:]
    sol_Ps=value.(Ps)[:]
    sol_Pg=value.(Pg)[:]
    sol_E=value.(E)[:,2]
    dPd_star = (Vector(sol_Pd).*ps.baseMVA)./Pd_max[:,1] # % load served
    dPs_star = Vector(sol_Ps).*ps.baseMVA
    dPg_star = Vector(sol_Pg).*ps.baseMVA
    dE_star = Vector(sol_E).*ps.baseMVA
    # add changes ps/psi structure
    ps.shunt.status = dPd_star;
    ps.storage.Ps = dPs_star;
    ps.storage.E = dE_star;
    ps.gen.Pg = dPg_star;
    @assert abs(sum(Pd_max[:,1].*ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg))<=2*tolerance
    @assert sum(ps.storage.E .< -tolerance)==0
    return ps
end

function crisp_mh_lsopf!(ps,dt,ug,ul,load_shed_cost;t_win=dt,w_g=0.1)
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
    Pdmax = (ps.shunt.P ./ ps.baseMVA)
    Pd1 = (ps.shunt.P ./ ps.baseMVA) .* ps.shunt.status
    if any(D.<1) || any(D.>n)
        error("Bad indices in shunt matrix")
    end
    # gen data
    ng = size(ps.gen.Pg,1)
    G = bi[ps.gen.bus]
    G_bus = sparse(G,collect(1:ng),1.,n,ng)
    Pg1 = (ps.gen.Pg ./ ps.baseMVA) .* ps.gen.status
    Pg_max = (ps.gen.Pmax ./ ps.baseMVA)
    RR = ps.gen.ramp30 ./30 .* dt ./ ps.baseMVA
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
    m = Model(with_optimizer(Gurobi.Optimizer))
    #m = Model(with_optimizer(Cbc.Optimizer))
    # variables
    @variable(m, Pd[1:nd]) # demand
    @variable(m, Pg[1:ng]) # generation
    @variable(m, Ps[1:ns]) # power flow into or out of storage (negative flow = charging)
    @variable(m, E[1:ns, 1:2]) # energy level in battery
    # fix battery starting charge
    for s in 1:ns
        fix(E[s,1], E1[s], force = true)
    end
    # variable bounds constraints
    @constraint(m, stPdcon, 0.0 .<= Pd[:] .<= Pdmax[:,1]) # load served limits
    @constraint(m, stPscon, Ps_min .<= Ps[:] .<= Ps_max) # storage power flow
    @constraint(m, stEPscon, E[:,2] .== (E[:,1] - ((dt/60) .* (Ps[:])))) # storage energy at next time step
    @constraint(m, stEcon, min_bat_E .<= (E[:,2]) .<= E_max) # storage energy
    @constraint(m, genPgucon, Pg[:] .<= ug[:,1] .* Pg_max[:,1]) # generator power limits upper
    @constraint(m, genPglcon, 0 .<= Pg[:]) # generator power limits lower
    if n > 1
        @variable(m, Theta[1:n])
        @constraint(m, Theta[1] .== 0); # set first bus as reference bus: V angle to 0
        brst = falses(nb);
        for b in 1:nb
            if ul[b,1] == 1
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
        @constraint(m, B*Theta[:] .== G_bus*Pg[:]+S_bus*Ps[:]-D_bus*Pd[:])
        # power flow limits
        @constraint(m, -flow_max[brst] .<= Xinv .* (Theta[F] - Theta[T]) .<= flow_max[brst])
    else
        @constraint(m, PBnoTheta, 0.0 .== G_bus*Pg[:]+S_bus*Ps[:]-D_bus*Pd[:])
    end
    # objective
    @objective(m, Max, sum(Pd) + sum(w_g.*Pg)) #TODO: constants
    ## SOLVE! ##
    optimize!(m)
    sol_Pd=value.(Pd)[:]
    sol_Ps=value.(Ps)[:]
    sol_Pg=value.(Pg)[:]
    sol_E=value.(E)[:,2]
    dPd_star = (Vector(sol_Pd).*ps.baseMVA)./ps.shunt.P # % load served
    dPs_star = Vector(sol_Ps).*ps.baseMVA
    dPg_star = Vector(sol_Pg).*ps.baseMVA
    dE_star = Vector(sol_E).*ps.baseMVA
    # add changes ps/psi structure
    ps.shunt.status = dPd_star;
    ps.storage.Ps = dPs_star;
    ps.storage.E = dE_star;
    ps.gen.Pg = dPg_star;
    @assert abs(sum(ps.shunt.P.*ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg))<=2*tolerance
    @assert sum(ps.storage.E .< -tolerance)==0
    return ps
end


function crisp_dcpf_NENY!(ps)
    # constants
    tolerance = 1e-4
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
    ng = size(ps.gen,1)
    G = bi[ps.gen.bus]
    Pg = ps.gen.P ./ ps.baseMVA .* ps.gen.status
    if any(G.<1) || any(G.>n)
        error("Bad indices in gen matrix")
    end
    Pg_bus = Array(sparse(G,ones(size(G)),Pg,n,1))
    # branch data
    brst = (ps.branch.status.==1)
    F = bi[ps.branch[brst,:f]]
    T = bi[ps.branch[brst,:t]]
    Xinv = (1 ./ ps.branch.X[brst])
    Bdc = sparse(F,T,-Xinv,n,n) + sparse(T,F,-Xinv,n,n) +
          sparse(T,T,+Xinv,n,n) + sparse(F,F,+Xinv,n,n)
    #find reference bus
    if isempty(ps.bus.id[ps.bus.kind.==3]) # note in Pavan's ps structure I beleive it's called 'kind' not 'bus_type'
        if isempty(ps.gen) || isempty(ps.shunt)
            ps.branch.Pf[brst] .= 0
            ps.branch.Pt[brst] .= 0
            ps.branch.Qf[brst] .= 0
            ps.branch.Qt[brst] .= 0
            if isempty(ps.shunt) && !isempty(ps.gen)
                ps.gen.P .= 0.0
            elseif !isempty(ps.shunt) && isempty(ps.gen)
                ps.shunt.P .= 0.0
            end
            return ps #TODO
        else
            maxGen = findmax(ps.gen.P)[2]
            busID = ps.gen[maxGen,:bus];
            isref = (busID.==ps.bus.id)
            nonref = .~isref
        end
    else
        isref = (ps.bus.kind.==3)
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
        is_refgen = (ps.gen.bus.==refbusid)
        if sum(is_refgen) != 1
            println("Multiple Gen on ref bus, splitting mismatch among them.")
        #if sum(is_refgen) != 1
        #    error("Must be exactly one ref generator")
        #end
            ps.gen.P[is_refgen] .-= (mismatch.*ps.baseMVA)/sum(is_refgen)
        else
            ps.gen.P[is_refgen] .-= (mismatch.*ps.baseMVA)
        end
    end
    # check the mismatch
    mis_check = sum(ps.gen.P.*ps.gen.status) - sum(ps.shunt.P.*ps.shunt.status)
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

function crisp_lsopf_bs!(ps,dt,ul,Pd_max,Pg_max1,load_shed_cost;t_win=dt,w_sl=10,w_g=0.1)
    timeline = 0:dt:t_win
    Ti = size(timeline,1)
    # constants
    tolerance = 1e-2
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    bi = ps.bi
    # load data
    nd = size(ps.shunt,1)
    D = bi[ps.shunt.bus]
    D_bus = sparse(D,collect(1:nd),1.,n,nd)
    Pdmax = (Pd_max[:,1] ./ ps.baseMVA)
    if any(D.<1) || any(D.>n)
        error("Bad indices in shunt matrix")
    end
    # gen data
    ng = size(ps.gen.Pg,1)
    Psl_max = zeros(ng);
    g_loads = (((ps.gen.state .== On) .| (ps.gen.state .== ShuttingDown) .| (ps.gen.state .== Damaged)
                .| (ps.gen.state .== WarmingUp)) .& .!ps.gen.black_start)
    Psl_max[g_loads] .= (ps.gen.service_load[g_loads] ./ ps.baseMVA)
    G = bi[ps.gen.bus]
    G_bus = sparse(G,collect(1:ng),1.,n,ng)
    Pg1 = (ps.gen.Pg ./ ps.baseMVA) .* ps.gen.status
    Pg_max = (Pg_max1 ./ ps.baseMVA) .* ps.gen.status
    ug = ps.gen.state .== On
    if sum(names(ps.gen) .== :RampRateMWMin) .!=0
        RR = ps.gen.RampRateMWMin .* dt ./ ps.baseMVA
    else
        RR = ps.gen.ramp30 ./30 .* dt ./ ps.baseMVA
    end
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
    m = Model(with_optimizer(Gurobi.Optimizer))
    #m = Model(with_optimizer(Cbc.Optimizer))
    # variables
    @variable(m, Pd[1:nd]) # demand
    @variable(m, Pg[1:ng]) # generation
    @variable(m, Psl[1:ng]) # generation
    @variable(m, Ps[1:ns]) # power flow into or out of storage (negative flow = charging)
    @variable(m, E[1:ns, 1:2]) # energy level in battery
    # fix battery starting charge
    for s in 1:ns
        fix(E[s,1], E1[s], force = true)
    end
    # variable bounds constraints
    @constraint(m, stPdcon, 0.0 .<= Pd[:] .<= Pdmax[:,1]) # load served limits
    @constraint(m, stPscon, Ps_min .<= Ps[:] .<= Ps_max) # storage power flow
    @constraint(m, stEPscon, E[:,2] .== (E[:,1] - ((dt/60) .* (Ps[:])))) # storage energy at next time step
    @constraint(m, stEcon, min_bat_E .<= (E[:,2]) .<= E_max) # storage energy
    @constraint(m, genPgucon, Pg[:] .<= ug[:] .* Pg_max) # generator power limits upper
    @constraint(m, genPssucon, Psl[:] .<= Psl_max) # generator power service load upper
    @constraint(m, genPglcon, 0 .<= Pg[:]) # generator power limits lower
    if n > 1
        @variable(m, Theta[1:n])
        @constraint(m, Theta[1] .== 0); # set first bus as reference bus: V angle to 0
        brst = falses(nb);
        for b in 1:nb
            if ul[b,1] == 1
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
        @constraint(m, B*Theta[:] .== G_bus*Pg[:]-G_bus*Psl[:]+S_bus*Ps[:]-D_bus*Pd[:])
        # power flow limits
        @constraint(m, -flow_max[brst] .<= Xinv .* (Theta[F] - Theta[T]) .<= flow_max[brst])
    else
        #one bus power balance
        @constraint(m, PBnoTheta, 0.0 .== G_bus*Pg[:]-G_bus*Psl[:]+S_bus*Ps[:]-D_bus*Pd[:])
    end
    # objective
    @objective(m, Max, w_sl*sum(Psl) + sum(Pd) + sum(w_g.*Pg)) #TODO: constants
    ## SOLVE! ##
    optimize!(m)
    sol_Pd=value.(Pd)[:]
    sol_Ps=value.(Ps)[:]
    sol_Pg=value.(Pg)[:]
    sol_Psl=value.(Psl)[:]
    sol_E=value.(E)[:,2]
    @assert abs(sum(sol_Pg)+sum(sol_Ps)-sum(sol_Pd)-sum(sol_Psl)) <= tolerance
    @assert abs(sum(sol_Pg.*ps.baseMVA)+sum(sol_Ps.*ps.baseMVA)-sum(sol_Pd.*ps.baseMVA)-sum(sol_Psl.*ps.baseMVA)) <= tolerance
    dPd_star = (Vector(sol_Pd).*ps.baseMVA)#./Pd_max
    dPs_star = Vector(sol_Ps).*ps.baseMVA
    dPg_star = Vector(sol_Pg).*ps.baseMVA
    dPsl_star = Vector(sol_Psl).*ps.baseMVA#./ps.gen.service_load
    dPsl_star[ps.gen.service_load.!=0] .= (dPsl_star[ps.gen.service_load.!=0]./ps.gen.service_load[ps.gen.service_load.!=0])
    dPsl_star[ps.gen.service_load.==0] .= 1# % service load served
    dE_star = Vector(sol_E).*ps.baseMVA
    @assert abs(sum(dPg_star)+sum(dPs_star)-sum(dPd_star)-sum(dPsl_star.*ps.gen.service_load)) <= tolerance
    # add changes ps/psi structure
    ps.shunt.status = dPd_star./Pd_max[:,1]; # % load served
    ps.storage.Ps = dPs_star;
    ps.storage.E = dE_star;
    ps.gen.Pg = dPg_star;
    ps.gen.sl_status = dPsl_star# % service load served
    @assert abs(sum(Pd_max[:,1].*ps.shunt.status)+sum(ps.gen.service_load.*ps.gen.sl_status)-sum(ps.storage.Ps)-sum(ps.gen.Pg))<=4*tolerance
    @assert sum(ps.storage.E .< -tolerance)==0
    return ps
end

#end
