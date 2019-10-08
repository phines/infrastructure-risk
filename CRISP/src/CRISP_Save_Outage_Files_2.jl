using CSV
include("CRISP_initiate.jl")
include("CRISP_network_gen.jl")
include("CRISP_LSOPF_gen1.jl")
#julia> rng = MersenneTwister(100);
#julia> Outages(1000,"data\\saved_ps\\case73_noPWS_lx2_n-1")
function Outages(Num,ps_folder;param_file = "",cascade=true)
    #constants
    debug=1;
    tolerance1 = 10^(-6);
    mult_factor = 1.5
    ## Num = number of failure scenarios to run through
    # initialize vector of costs from events
    ## load the case data
    ps = import_ps("$ps_folder")
    ps.shunt = ps.shunt[ps.shunt.P .!=0.0,:]
    crisp_dcpf_g_s!(ps)
    if isempty(param_file)
        # parameters of distributions for line outages and recovery times
        s_line = 2.56;#lines_dist[1];
        maxLinesOut = length(ps.branch.f); # => k in zipf distribution
        mu_line = 3.66;#lines_dist[2];
        sigma_line = 2.43;#lines_dist[3];
        lambda_gen = 1;
    end
    for iterat in 1:Num
        # step 1
        if cascade
            Lines_Init_State = line_state_cascade!(ps,s_line,maxLinesOut,mu_line,sigma_line)
            dt = 10
            subgraph = find_subgraphs(ps);
            M = Int64(findmax(subgraph)[1]);
            ps_islands = build_islands(subgraph,ps);
            for i in 1:M
                psi = ps_subset(ps,ps_islands[i]);
                # run lsopf
                crisp_lsopf_g_s!(psi,dt);
                ps.gen[ps_islands[i].gen,:Pg] = psi.gen.Pg
                ps.shunt[ps_islands[i].shunt,:status] = psi.shunt.status
                ps.storage[ps_islands[i].storage,:E] = psi.storage.E
                ps.storage[ps_islands[i].storage,:Ps] = psi.storage.Ps
            end
            line_rec_times_add_comms!(ps,Lines_Init_State,mult_factor)
        else
            Lines_Init_State = line_state!(ps,s_line,maxLinesOut,mu_line,sigma_line)
        end
        Gens_Init_State = gen_state!(ps,lambda_gen,mu_line,sigma_line)
        if debug==1
            CSV.write("data\\outage_data\\communication_factor\\out_case73_noPWS_lx2_n-1_lines$iterat.csv", Lines_Init_State)
            CSV.write("data\\outage_data\\communication_factor\\out_case73_noPWS_lx2_n-1_gens$iterat.csv", Gens_Init_State)
        end
    end
end
#=
#for debugging
include("src\\CRISP_initiate.jl")
include("src\\CRISP_network_gen.jl")
include("src\\CRISP_LSOPF_gen1.jl")
=#
function line_rec_times_add_comms!(ps,Lines_Init_State,mult_factor)
    line_status  = Lines_Init_State[:,1]
    r_time = Lines_Init_State[:,2]
    nl = size(ps.branch.f)[1]
    r_time_new = zeros(nl)
    for b in 1:nl
        if line_status[b] .!= 1
            fr = ps.branch.f[b]
            to = ps.branch.t[b]
            if (sum(fr .== ps.shunt.bus) != 0) || (sum(to .== ps.shunt.bus) != 0)
                if sum(ps.shunt.status[ps.shunt.bus .== fr] .!= 1) != 0
                    r_time_new[b] = r_time[b].*mult_factor
                    println("comms effect recovery time of line $b")
                elseif sum(ps.shunt.status[ps.shunt.bus .== to] .!= 1) != 0
                    r_time_new[b] = r_time[b].*mult_factor
                    println("comms effect recovery time of line $b")
                else
                    r_time_new[b] = r_time[b]
                end
            else
                r_time_new[b] = r_time[b]
            end
        end
    end
    Lines_Init_State = DataFrame(state = line_status, recovery_time = r_time_new)
    return ps, Lines_Init_State
end
