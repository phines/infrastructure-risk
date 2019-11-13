using SpecialFunctions
using Random
include("CRISP_network_gen.jl")
# interaction functions

function communication_interactions(ps,restoration_times,comm_count,t;n_sigmoids = 5,param=[1 2 3 4 5; 4 16 24 48 48])
    if comm_count > n_sigmoids
    else
        f = sigmoid(param[:,comm_count],t)
        if rand(rng,1) < f
            comm_count += 1
        end
    end
    return restoration_times
end

function logistic(param,t;L=1)
    f = L/(1+exp(-param[1]*(t-param[2])))
    return f
end
