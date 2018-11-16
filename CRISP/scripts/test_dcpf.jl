include("../src/CRISP_LSOPF.jl")
include("../src/parser.jl")

# load the case data
ps = mp2ps("../data/case6ww.m")

# dcpf
crisp_dcpf!(ps)
