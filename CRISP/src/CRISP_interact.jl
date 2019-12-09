using SpecialFunctions
using Random
include("CRISP_network.jl")
# interaction functions
#constants
four_days = 4*24
week = 7*24
two_weeks = 14*24;
three_months = 90*24;
function natural_gas_interactions!(ps,Lines_Init_State,Gens_Init_State;
                                range_a=two_weeks,range_b=three_months)
    if rand(rng,1)[1] >= sum(Lines_Init_State.state)/length(Lines_Init_State.state)
        println("ng_fails = 1")
        ng_rest_time = rand(rng,collect(range_a:range_b),1)
        ng = size(ps.gen,1)
        state = Gens_Init_State.state
        recovery_times = Gens_Init_State.recovery_time
        gen = collect(1:ng)
        ng_gen = gen[ps.gen.Fuel .== "Natural Gas"]
        for g in ng_gen
            state[g] .= 0
            recovery_times[g] .= ng_rest_time
        end
        return Gens_Init_NG_State = DataFrame(state = state, recovery_time = recovery_times)
    else
        return Gens_Init_State
    end
end

function nuclear_poissoning(ps,n; time_range = four_days:week)
    r = rand(rng,tim_range,n)
    return r
end

function compound_rest_times!(ps,Pdmax,restoration_times,factor,ti)
    tolerance = 10^(-8)
    println(ti)
    loadshed = ps.shunt.status.*Pdmax
    prob_check = loadshed .>= rand(rng,length(loadshed))
    buses_effected = ps.shunt.bus[prob_check]
    for b in buses_effected
        lines_f = ps.branch.f .== b
        restoration_times[lines_f] = restoration_times[lines_f].*factor
        lines_t = ps.branch.t .== b
        restoration_times[lines_t] = restoration_times[lines_t].*factor
    end
    return restoration_times
end

function communication_interactions(ps,restoration_times,comm_battery_limits,t,factor)
    l = restoration_times .> 0
    B = Int64.(ps.bus.id)
    F = Int64.(ps.branch.f[l])
    T = Int64.(ps.branch.t[l])
    CBLF = Int64.(zeros(length(F)))
    CBLT = Int64.(zeros(length(T)))
    LSF = Int64.(zeros(length(F)))
    LST = Int64.(zeros(length(T)))
    for i in 1:length(F)
        CBLF[i] = comm_battery_limits[F[i] .== B][1]
        CBLT[i] = comm_battery_limits[T[i] .== B][1]
        if sum(ps.shunt.bus .== F[i]) >= 1
            LSF[i] += ps.shunt.status[ps.shunt.bus .== F[i]][1]
        end
        if sum(ps.shunt.bus .== T[i]) >= 1
            LST[i] += ps.shunt.status[ps.shunt.bus .== T[i]][1]
        end
        if CBLF[i] == t
            g = rand(rng,1) < LSF
            if g
                restoraion_times[i] = restoraion_times[i].*factor
            end
        elseif CBLT[i] == t
            if rand(rng,1) < LST
                restoraion_times[i] = restoraion_times[i].*factor
            end
        end
    end
    return restoration_times
end

function comm_battery_lim(n,a,b)
    #find the battery limits of comms, a and b are the range, and n are the number of buses.
    return comm_bl = rand(rng,a:b,n);
end

function communication_logistic_interactions(ps,restoration_times,comm_count,t;
                    n_sigmoids = 5,param=[0.1 0.2 0.3 0.4 0.5; 4 16 24 48 72])
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
