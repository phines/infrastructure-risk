#set up packages
using CSV; using DataFrames; using SpecialFunctions;
include("CRISP_LSOPF_1.jl")
function RLSOPF!(ps,failures,recovery_times)
    times = sort(recovery_times)
    load_shed = zeros(length(times))
    # set load shed for the step just before restoration process
    load_shed[1] = sum(ps.shunt[:P]);
    for i = 1:length(times)
        T = times[i];
        # set failed branches to status 0
        failures[T.==recovery_times] = 0;
        # apply to network
        ps.branch[:,:status] = failures;
        # run the dcpf
        crisp_dcpf!(ps)
        # run lsopf
        (dPd, dPg) = crisp_lsopf(ps)
        # apply the results
        ps.gen[:Pg]  += dPg
        ps.shunt[:P] += dPd
        crisp_dcpf!(ps)
        # set load shed for this time step
        load_shed[i+1] = sum(ps.shunt[:P]);
    end
    Restore = DataFrame(time = times, load_shed = load_shed)
    return Restore
end
