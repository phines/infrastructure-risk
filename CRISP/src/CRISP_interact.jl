using SpecialFunctions
using Random
include("CRISP_network_gen.jl")
# interaction functions

function communication_interactions(ps,restoration_times,comm_battery_limits,t)
    l = restoration_times > 0;
    B = ps.bus.id;
    F = ps.branch.f[l];
    T = ps.branch.t[l];
    BL = falses(length(comm_battery_limits));
    for i in 1:length(F)
        BL[F[i] .== B] .= true;
        BL[T[i] .== B] .= true;
    end

    if sum(comm_battery_limits[BL] .== t) > 0
        rand(rng,1) <= load_shed
    end
    return restoration_times
end

function communication_logistic_interactions(ps,restoration_times,comm_count,t;
                        n_sigmoids = 5,param=[1 2 3 4 5; 4 16 24 48 72])
    if comm_count > n_sigmoids
    else
        f = logistic(param[:,comm_count],t)
        if rand(rng,1) < f
            comm_count += 1
        end
    end
    return restoration_times, comm_count
end

function logistic(param,t;L=1)
    f = L/(1+exp(-param[1]*(t-param[2])))
    return f
end

function natural_gas_interactions(ps)
    return ps
end

function nuclear_poissoning(ps)
    return ps
end
