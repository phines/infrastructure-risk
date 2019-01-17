module CRISP
#set up packages
using CSV; using DataFrames; using SpecialFunctions;
include("s1-initiate1.jl")
include("CRISP_LSOPF_1.jl")
include("CRISP_RLSOPF.jl")

end
