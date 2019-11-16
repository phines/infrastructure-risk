using SpecialFunctions
using Random
include("CRISP_network_gen.jl")
# interaction functions

function communication_interactions(ps,restoration_times,comm_battery_limits,t,factor)
    l = restoration_times > 0;
    B = ps.bus.id;
    F = ps.branch.f[l];
    T = ps.branch.t[l];
    CBLF = zeros(length(F))
    CBLT = zeros(length(T))
    LSF = zeros(length(F))
    LST = zeros(length(T))
    for i in 1:length(F)
        CBLF[i] = comm_battery_limits[F[i] .== B];
        CBLT[i] = comm_battery_limits[T[i] .== B];
        LSF[i] += ps.load.status[ps.load.bus .== F[i]]
        LST[i] += ps.load.status[ps.load.bus .== T[i]]
        if CBLF[i] == t
            g = rand(rng,1) < LSF
            if g
                restoraion_times[i] = restoraion_times[i].*factor
            end
        end
        if CBLT[i] == t & !g
            if rand(rng,1) < LST
                restoraion_times[i] = restoraion_times[i].*factor
            end
        end
    end
    return restoration_times
end
function comm_battery_limits (n,a,b)
    return comm_bl = rand(rng,a:b,n);
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

function nuclear_poissoning(ps,n; time_range = (4*24):(7*24))
    r = rand(rng,tim_range,n)
    return r
end
