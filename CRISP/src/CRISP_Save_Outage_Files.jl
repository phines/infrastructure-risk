using CSV
include("CRISP_initiate.jl")
include("CRISP_network_gen.jl")
#julia> rng = MersenneTwister(100);
#julia> Outages(1000,"data\\saved_ps\\case73_noPWS_lx2_n-1")
function Outages(Num,ps_folder;param_file = "",cascade=false)
    debug=1;
    tolerance1 = 10^(-6);
    ## Num = number of failure scenarios to run through
    # initialize vector of costs from events
    ## load the case data
    ps = import_ps("$ps_folder")
    ps.shunt = ps.shunt[ps.shunt.P .!=0.0,:]
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
        else
            Lines_Init_State = line_state!(ps,s_line,maxLinesOut,mu_line,sigma_line)
        end
        Gens_Init_State = gen_state!(ps,lambda_gen,mu_line,sigma_line)
        if debug==1
            CSV.write("data\\outage_data\\1\\out_case73_noPWS_lx2_n-1_lines$iterat.csv", Lines_Init_State)
            CSV.write("data\\outage_data\\1\\out_case73_noPWS_lx2_n-1_gens$iterat.csv", Gens_Init_State)
        end
    end
end
#=
#for debugging
include("src\\CRISP_initiate.jl")
include("src\\CRISP_network_gen.jl")
=#
