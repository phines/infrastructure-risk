using CSV
include("..\\src\\CRISP_initiate.jl")
include("..\\src\\CRISP_LSOPF_gen_stor.jl")
include("..\\src\\CRISP_RLOPF_gen_stor.jl")
include("..\\src\\CRISP_RT.jl")
include("..\\src\\CRISP_network_gen.jl")

    ## load the case data
    #case = "data\\saved_ps\\case73_noPWS\\"
    case = "data\\saved_ps\\case39_n-1_gen\\"
    ps = import_ps("$case")
    ps.shunt = ps.shunt[ps.shunt.P .!=0,:];
    #add columns to keep track of the time each generator is on or off
    if sum(names(ps.gen).==:time_on) == 0
        ps.gen.time_on = zeros(length(ps.gen.Pg));
    end
    if sum(names(ps.gen).==:time_off) == 0
        ps.gen.time_off = zeros(length(ps.gen.Pg));
    end
    if sum(names(ps.gen).==:off) == 0
        ps.gen.off = zeros(length(ps.gen.Pg));
        ps.gen.off[ps.gen.Pg .== 0] .= 1;
    end
    crisp_opf_initiate!(ps,0.1)
    total = sum(ps.shunt[:P]);
    Pd_max = deepcopy(ps.shunt[:P]);
    ps0 = deepcopy(ps);
    Lines_Init_State = CSV.read("results\\experiments_gen_stor\\case39\\res_out_case39_n-1_gen IC7 lines.csv", allowmissing=:none)
    Gens_Init_State = CSV.read("results\\experiments_gen_stor\\case39\\res_out_case39_n-1_gen IC7 gens.csv", allowmissing=:none)
    #Lines_Init_State = CSV.read("results\\experiments_gen_stor\\res_out_case73_noPWS IC9 lines.csv", allowmissing=:none)
    #Gens_Init_State = CSV.read("results\\experiments_gen_stor\\res_out_case73_noPWS IC9 gens.csv", allowmissing=:none)
        l_failures = Lines_Init_State[:,1];
        ps.branch.status[l_failures .==0] .= 0;
        NumLinesOut = length(l_failures) - sum(l_failures);
        l_recovery_times = Lines_Init_State[:,2];
        # generator states and recovery times
        g_failures = Gens_Init_State[:,1];
        ps.gen.status[g_failures.==0] .= 0;
        ps.gen.Pg[g_failures.==0] .= 0;
        g_recovery_times = Gens_Init_State[:,2];
        #check for islands
        subgraph = find_subgraphs(ps);
        M = Int64(findmax(subgraph)[1]);
        ps_islands = build_islands(subgraph,ps);
        for i in 1:M
            psi = ps_subset(ps,ps_islands[i]);
            dt = 1;
            # run lsopf
            crisp_lsopf_g_s!(psi,dt);
            @assert 10^(-6)>=abs(sum(psi.shunt.P .*psi.shunt.status)-sum(psi.storage.Ps)-sum(psi.gen.Pg))
            ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
            ps.shunt[ps_islands[i].shunt,:status] = psi.shunt.status
            ps.storage[ps_islands[i].storage,:E] = psi.storage.E
            ps.storage[ps_islands[i].storage,:Ps] = psi.storage.Ps
        end
        @assert 10^(-6)>=abs(sum(ps.shunt.P .*ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg))
        @assert total>=sum(ps.shunt.P)
        @assert total>=sum(ps.shunt.P .*ps.shunt.status)
        @assert 10^(-6)>=abs(sum(ps.shunt.P .*ps.shunt.status)-sum(ps.storage.Ps)-sum(ps.gen.Pg))
        total-sum(ps.shunt.P .*ps.shunt.status)
        ## run step 3
        dt = 1
        ti = dt;#10
        t0 = 10
        #crisp_mh_rlopf!(ps,dt,time)
        Restore = crisp_Restore(ps,l_recovery_times,g_recovery_times,dt,ti,t0)
        ###find the time to restore the grid to 99.9% load served
        K = abs.(Restore.perc_load_served .- 1) .<= 0.001;
        K[1] = false;
        if isempty( Restore.time[K]) && (size(Restore)[1] > 1)
            error("Never solved to full restoration.")
        elseif isempty( Restore.time[K]) && (size(Restore)[1] == 1)
            LoadServedTime = t0;
        else
            R = Restore.time[K];
            LoadServedTime = R[1] - Restore.time[1];
            ## run step 4
            ResilienceTri = crisp_res(Restore);
        end


#=
#for debugging
include("src\\CRISP_initiate.jl")
include("src\\CRISP_LSOPF_gen_stor.jl")
include("src\\CRISP_RLOPF_gen_stor.jl")
include("src\\CRISP_RT.jl")
include("src\\CRISP_network_gen.jl")
=#
