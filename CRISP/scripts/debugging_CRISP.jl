# Debugging results/case39/*case39_9.*

using CSV

include("..\\src\\CRISP_LSOPF.jl")
include("..\\src\\CRISP_RLSOPF.jl")
include("..\\src\\CRISP_RT.jl")
include("..\\src\\CRISP_network_segments.jl")
include("..\\src\\parser3.jl")
include("..\\src\\export_ps.jl")
#include("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\src\\parser.jl")
## load the case data
ps = mp2ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case39.m")
crisp_dcpf!(ps)
total = sum(ps.shunt[:P]);
Pd_max = deepcopy(ps.shunt[:P]);

#save restoration picture to:
filename = "Fixed_case39_9";

#pull failures from csv file
Lines_Init_State = CSV.read("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\results\\case39\\test_initial_outage_case39_9.csv");
state = Lines_Init_State[:,1];
recovery_times = Array{Float64}(Lines_Init_State[1:end,2]);
failures = Array{Int64}(state[1:end]);
ps.branch.status = failures;

#check for islands
subgraph = find_subgraphs(ps);
M = Int64(findmax(subgraph)[1]);
ps_islands = build_islands(subgraph,ps);
for i in 1:M
    psi = ps_subset(ps,ps_islands[i]);
    n = size(psi.bus,1) # the number of buses
    if n==1
        if (!isempty(psi.gen) && isempty(psi.shunt)) || (isempty(psi.gen) && !isempty(psi.shunt))
            dPg = -(psi.gen[:Pg] ./ psi.baseMVA .* psi.gen[:status]);
            dPd = -(psi.shunt[:P] ./ psi.baseMVA .* psi.shunt[:status]);
            dPd_star = dPd.*ps.baseMVA;
            dPg_star = dPg.*ps.baseMVA;
            ps.gen[ps_islands[i].gen,:Pg]  += dPg_star;
            ps.shunt[ps_islands[i].shunt,:P] += dPd_star;
        elseif !isempty(psi.gen) && !isempty(psi.shunt)
            Pd = psi.shunt[:P] ./ psi.baseMVA .* psi.shunt[:status]
            Pg = psi.gen[:Pg] ./ psi.baseMVA .* psi.gen[:status]
            Pg_cap = psi.gen[:Pmax] ./ psi.baseMVA .* psi.gen[:status]
            if Pg_cap >= Pd
                dPd = 0;
                dPg = Pd-Pg;
            else
                dPd = Pd-Pg_cap;
                dPg = Pg_cap-Pg;
            end
            dPd_star = dPd.*ps.baseMVA;
            dPg_star = dPg.*ps.baseMVA;
            ps.gen[ps_islands[i].gen,:Pg]  = dPg_star;
            ps.shunt[ps_islands[i].shunt,:P] = dPd_star;
        else
            ps.gen[ps_islands[i].gen,:Pg]  = psi.shunt[:P].*0
            ps.shunt[ps_islands[i].shunt,:P] = psi.gen[:Pg].*0
        end
    else
        ## run step 2
        # run the dcpf
        crisp_dcpf!(psi);
        # run lsopf
        (dPd, dPg) = crisp_lsopf(psi);
        # apply the results
        ps.gen[ps_islands[i].gen,:Pg]  += dPg;
        ps.shunt[ps_islands[i].shunt,:P] += dPd;
        crisp_dcpf!(psi);
    end
end

# save ps structure
export_ps(ps,"data\\saved_ps\\case39_severe")

## run step 3
Restore = RLSOPF!(total,ps,failures,recovery_times,Pd_max);#,load_cost) # data frame [times, load shed in cost per hour]
#make_plots(Restore,filename)
using Plots; using StatsPlots
plot1 = @df Restore plot(:time, :perc_load_served,
      title = "Resilience Triangle",
      xlabel = "time", ylabel = "load served (%)")
plot2 = @df Restore plot(:time, :num_lines_out,
            title = "Line Restoration",
            xlabel = "time", ylabel = "number of lines out")
#putting 2 plots together
P = plot(plot1,plot2,layout = (1,2))
# save a png
png(P,"results\\$filename")
## run step 4
ResilienceTri = crisp_res(Restore);
