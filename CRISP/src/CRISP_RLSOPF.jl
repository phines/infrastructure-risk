#set up packages
using CSV; using DataFrames; using SpecialFunctions;
include("CRISP_LSOPF_1.jl")
function RLSOPF!(ps,failures,recovery_times)
    times = sort(recovery_times)
    for i = 1:length(times)
        
    end
end
