using SpecialFunctions
using Random
using Distributions
include("CRISP_network.jl")

# interaction functions
#constants
four_days = 4*24*60
week = 7*24*60
two_weeks = 14*24*60;
three_months = 90*24*60;

function natural_gas_interactions!(ps,Lines_Init_State,Gens_Init_State;
                                range_a=two_weeks,range_b=three_months)
    if rand(rng,1)[1] >= sum(Lines_Init_State.state)/length(Lines_Init_State.state)
        println("ng_fails = 1")
        ng_rest_time = rand(rng,collect(range_a:range_b),1)
        ng = size(ps.gen,1)
        state = Gens_Init_State.state
        recovery_times = Gens_Init_State.recovery_time
        gen = collect(1:ng)
        ng_gen = gen[ps.gen.Fuel .== "NG"]
        for g in ng_gen
            state[g] = 0
            recovery_times[g] = ng_rest_time[1]
        end
        return Gens_Init_NG_State = DataFrame(state = state, recovery_time = recovery_times)
    else
        return Gens_Init_State
    end
end

function natural_gas_int_lognorm!(ps,Lines_Init_State,Gens_Init_State;
                                mu = 3, sigma = 1)
    lognorm = LogNormal(mu,sigma) # in days
    if rand(rng,1)[1] >= sum(Lines_Init_State.state)/length(Lines_Init_State.state)
        println("ng_fails = 1")
        ng_rest_time = rand(rng,lognorm,1) # in days
        ng = size(ps.gen,1)
        state = Gens_Init_State.state
        recovery_times = Gens_Init_State.recovery_time
        gen = collect(1:ng)
        ng_gen = gen[ps.gen.Fuel .== "NG"]
        for g in ng_gen
            state[g] = 0
            recovery_times[g] = ng_rest_time[1].*60*24 # from days to minutes
        end
        return Gens_Init_NG_State = DataFrame(state = state, recovery_time = recovery_times)
    else
        return Gens_Init_State
    end
end

function nuclear_pois_lognorm!(ps,Pg_i,g_recovery_times,ti; mu = 1, sigma = 0.5)
    lognorm = LogNormal(mu,sigma) # days
    tolerance = 10^(-4)
    ng = size(ps.gen,1)
    gen = collect(1:ng)
    nuc_gen = gen[ps.gen.Fuel .== "Nuclear"]
    for g in nuc_gen
        if (Pg_i[g,1] > ps.gen.Pg[g]) .& (Pg_i[g,1] > Pg_i[g,2]) .& (Pg_i[g,1] > Pg_i[g,3]) .&(abs(ps.gen.Pg[g]) <= tolerance)
            ps.gen.status[g] = 0
            if g_recovery_times[g] .== 0
                g_recovery_times[g] = ti + (rand(rng,lognorm,1)[1].*60*24)  # into units of minutes
            else
                g_recovery_times[g] += (rand(rng,lognorm,1)[1].*60*24 )
            end
        end
    end
    return g_recovery_times
end

function nuclear_poissoning!(ps,Pg_i,g_recovery_times,ti; time_range = four_days:week)
    tolerance = 10^(-4)
    ng = size(ps.gen,1)
    gen = collect(1:ng)
    nuc_gen = gen[ps.gen.Fuel .== "Nuclear"]
    for g in nuc_gen
        if (Pg_i[g] > ps.gen.Pg[g]) .& (abs(ps.gen.Pg[g]) <= tolerance)
            ps.gen.status[g] = 0
            if g_recovery_times[g] .== 0
                g_recovery_times[g] = ti+rand(rng,time_range,1)[1]
            else
                g_recovery_times[g] += rand(rng,time_range,1)[1]
            end
        end
    end
    return g_recovery_times
end

function compound_rest_times!(ps,restoration_times,factor,ti)
    tolerance = 10^(-8)
    prob_check = ps.shunt.status .< 1#rand(rng,length(ps.shunt.status))
    buses_effected = ps.shunt.bus[prob_check]
    for b in buses_effected
        lines_f = (ps.branch.f .== b)
        restoration_times[lines_f] .*= factor
        lines_t = (ps.branch.t .== b)
        restoration_times[lines_t] .*= factor
    end
    return restoration_times
end

function communication_interactions!(ps,restoration_times,comm_battery_limits,t,factor)
    println("INTO COMMMS INTERACTIONS")
    l = restoration_times .> 0
    r = restoration_times[l]
    B = Int64.(ps.bus.id)
    F = Int64.(ps.branch.f[l])
    T = Int64.(ps.branch.t[l])
    CBLF = Int64.(zeros(length(F)))
    CBLT = Int64.(zeros(length(T)))
    LSF = ones(length(F))
    LST = ones(length(T))
    for i in 1:length(F)
        CBLF[i] = comm_battery_limits[F[i] .== B][1]
        CBLT[i] = comm_battery_limits[T[i] .== B][1]
        if sum(ps.shunt.bus .== F[i]) >= 1
            LSF[i] -= ps.shunt.status[ps.shunt.bus .== F[i]][1]
            println("MADE IT TO CHECK THE LOAD SHED")
        end
        if sum(ps.shunt.bus .== T[i]) >= 1
            LST[i] -= ps.shunt.status[ps.shunt.bus .== T[i]][1]
            println("MADE IT TO CHECK THE LOAD SHED")
        end
        if CBLF[i] == (t/60)
            println("MADE IT TO CHECK THE BATTERY LIFE")
            if rand(rng,1)[1] > LSF[i]
                println(r[i])
                r[i] *= factor
                println(r[i])
                println("CHANGED THE RESTORATION TIME??????")
            end
        elseif CBLT[i] == (t/60)
            println("MADE IT TO CHECK THE BATTERY LIFE")
            if rand(rng,1)[1] > LST[i]
                println(r[i])
                println(factor)
                r[i] *= factor
                println(r[i])
                println("CHANGED THE RESTORATION TIME??????")
            end
        end
    end
    restoration_times[l] .= r
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
