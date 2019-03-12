include("../src/CRISP_LSOPF.jl")
#include("../src/CRISP_LSOPF_1.jl")
include("../src/parser.jl")

# load the case data
ps = ("../data/case6ww.m")
#ps = mp2ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case39.m")

# dcpf
crisp_dcpf!(ps)
