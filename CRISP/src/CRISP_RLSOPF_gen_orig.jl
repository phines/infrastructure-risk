using DataFrames;
using SpecialFunctions;
using JuMP;
using Clp;
using Gurobi;
using Cbc;
include("CRISP_LSOPF.jl")
include("CRISP_network.jl")

function RLSOPF_g!(ps,l_failures,g_failures,l_recovery_times,g_recovery_times,gen_startup,Pd_max;t0 = 10, load_cost=0)
    #add columns to keep track of the time each generator is on or off
    if sum(names(ps.gen).==:time_on) == 0
        ps.gen.time_on = zeros(length(ps.gen.Pg));
    end
    if sum(names(ps.gen).==:time_off) == 0
        ps.gen.time_off = zeros(length(ps.gen.Pg));
    end
    # constants
    deltaT = 5; # time step in minutes;
    tolerance = 1e-6
    if sum(load_cost)==0
        load_cost = ones(length(ps.shunt.P));
    end
    l_rec_t = l_recovery_times[l_recovery_times.!=0];
    g_rec_t = g_recovery_times[g_failures.==0] + gen_startup[g_failures.==0];
    g_rec = zeros(length(g_failures));
    g_rec[g_failures.==0] = g_rec_t;
    Time = sort([l_rec_t; g_rec_t])
    total_i = length(Time);
    load_shed = zeros(total_i+2);
    load_shed[1] = 0;
    lines_out = zeros(total_i+2);#zeros(length(l_times)+2);
    lines_out[2] = length(l_failures) - sum(l_failures);
    gens_out = zeros(total_i+2);#zeros(length(l_times)+2);
    gens_out[2] = length(g_failures) - sum(g_failures);
    # set load shed for the step just before restoration process
    load_shed[2] = sum(load_cost.*(Pd_max - ps.shunt[:P]));
    time_gen_off = ones(length(ps.gen.Pg));
    time_gen_off[ps.gen.Pg .!=0] .=0;
    for i = 1:length(Time)
        T = Time[i];
        # set failed branches to status 0
        l_failures[T.>=l_recovery_times] .= 1;
        lines_out[i+2] = length(l_failures) - sum(l_failures);
        # apply to network
        ps.branch[:,:status] = l_failures;
        # set newly available generators to status 1;
        current_gen_out = ps.gen.status;
        # update generators operational
        g_failures[T.>=g_rec] .= 1;
        ps.gen.status = g_failures;
        gens_out[i+2] = length(g_failures) - sum(g_failures);
        #check for islands
        subgraph = find_subgraphs(ps);# add Int64 here hide info here
        M = Int64(findmax(subgraph)[1]);
        if M==1
            crisp_dcpf_g!(ps);
            # run lsopf
            crisp_rlopf_g!(ps,Pd_max);
            crisp_dcpf_g!(ps);
        else
            ps_islands = build_islands(subgraph,ps);# at some point check for changes in islands and don't run power flows if no change
            ## for every island that changed (eventually)

            for j in 1:M
                psi = ps_subset(ps,ps_islands[j]);
                    # run the dcpf
                    crisp_dcpf_g!(psi);
                    # run lsopf
                    crisp_rlopf_g!(psi,Pd_max[ps_islands[j].shunt]);
                    # apply the results
                    ps.gen[ps_islands[j].gen,:Pg]  = psi.gen.Pg;
                    ps.shunt[ps_islands[j].shunt,:P] = psi.shunt.P;
                    crisp_dcpf_g!(psi);
            end
        end
        @assert abs(sum(ps.gen.Pg)-sum(ps.shunt.P))<=tolerance
        # set load shed for this time step
        load_shed[i+2] = sum(load_cost.*(Pd_max - ps.shunt[:P]));
    end
    times = [0.0;t0*1.0;Time.+t0*1.0];
    perc_load_served = (sum(load_cost.*Pd_max) .- load_shed)./sum(load_cost.*Pd_max);
    Restore = DataFrame(time = times, load_shed = load_shed, perc_load_served = perc_load_served, num_lines_out = lines_out, num_gens_out = gens_out);
    return Restore
end

function crisp_rlopf_g!(ps,Pd_max)
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
        gst = (ps.gen.status .== 1);
        ng = size(ps.gen[gst,:Pg],1)
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
        m1 = Model(with_optimizer(Gurobi.Optimizer))
        # variables
        @variable(m1,dPd[1:nd])
        @variable(m1,ndPg[1:ng])
        @variable(m1,pdPg[1:ng])
        @variable(m1,ug[1:ng],Bin)
        @variable(m1,dTheta[1:n])
        # variable bounds
        @constraint(m1,-Pd.<=dPd.<=(Pd_max./ps.baseMVA - Pd));
        @constraint(m1, ug.*(Pg_min) .<= Pg+ndPg)
        @constraint(m1, ug.*(Pg_max) .>= Pg+pdPg)
        @constraint(m1, pdPg .>= 0)
        @constraint(m1, ndPg .<= 0)
        @constraint(m1,dTheta[1] == 0)
        # objective
        @objective(m1,Max,sum(dPd)+0.01*(0.1*sum(ndPg)-0.1*sum(pdPg)+sum(ug))) # serve as much load as possible
        # mapping matrix to map loads/gens to buses
        # Power balance equality constraint
        @constraint(m1,B*dTheta .== G_bus*ndPg+G_bus*pdPg-D_bus*dPd);#M_G*dPg - M_D*dPd)
        # Power flow constraints
        @constraint(m1,-flow_max .<= flow0 + Xinv.*(dTheta[F] - dTheta[T]) .<= flow_max)
        ### solve the model ###
        optimize!(m1);
        # collect/return the outputs
        sol_dPd=value.(dPd)
        sol_ndPg=value.(ndPg)
        sol_pdPg=value.(pdPg)
        dPd_star = sol_dPd.*ps.baseMVA
        dPg_star = (sol_pdPg+sol_ndPg).*ps.baseMVA
        ps.shunt.P += dPd_star; #changes ps structure
        ps.gen.Pg[gst] += dPg_star; #changes ps structure
    else
        if ((!isempty(ps.gen.status.==1) && isempty(ps.shunt)) || (isempty(ps.gen.status.==1) && !isempty(ps.shunt))
             || (isempty(ps.gen.status.==1) && isempty(ps.shunt)))
            ps.gen.Pg  .= ps.gen.Pg.*0.0;
            ps.shunt.P .= ps.shunt.P.*0.0;
        elseif !isempty(ps.gen.status.==1) && !isempty(ps.shunt)
            gst = (ps.gen.status .== 1)
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
            ps.gen.Pg  = ps.gen.Pg.*0.0;
        end
    end
    #adding criteria that should produce errors if incorrect.
    @assert abs(sum(ps.shunt.P)-sum(ps.gen.Pg))<=2*tolerance
    #@assert 0.<=ps.shunt.P
    return ps
end
