include("../src/CRISP_LSOPF.jl")
include("../src/CRISP_network.jl")

ps = import_ps("../../../NE_NY_data")
crisp_dcpf_NENY!(ps)
