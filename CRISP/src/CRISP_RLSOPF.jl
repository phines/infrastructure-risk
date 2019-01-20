#set up packages
using CSV; using DataFrames; using SpecialFunctions;
include("CRISP_LSOPF_1.jl")
function RLSOPF!(totalp,ps,failures,recovery_times)
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
        (dPd, dPg) = crisp_lsopf(ps)
        # apply the results
        ps.gen[:Pg]  += dPg
        ps.shunt[:P] += dPd
        crisp_dcpf!(ps)
        # set load shed for this time step
        load_shed[i+1] = totalp - sum(ps.shunt[:P]);
    end
    times = [0;times];
    Restore = DataFrame(time = times, load_shed = load_shed)
    return Restore
end
