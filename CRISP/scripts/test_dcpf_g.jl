include("../src/CRISP_LSOPF_gen.jl")
include("../src/CRISP_network.jl")

# load the case data
#ps = import_ps("../data/case6ww/")
ps = import_ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case6ww\\")

# dcpf
crisp_dcpf_g!(ps)
