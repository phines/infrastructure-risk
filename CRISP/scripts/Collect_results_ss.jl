using CSV
using DataFrames
using Glob
include("../src/CRISP_initiate.jl")
include("../src/CRISP_network.jl")
case = "data\\saved_ps\\case73_noPWS_lx2_n-1"
ps = import_ps("..\\"*case)
TotalL = length(ps.branch.f)
TotalG = length(ps.gen.bus)
events = "..\\data\\strat_sample_bins\\casc_73bus"
Line = CSV.File(events*"/line_bins.csv") |> DataFrame
Gen = CSV.File(events*"/gen_bins.csv") |> DataFrame
L = length(Line.line_bins)
G = length(Gen.gen_bins)
Line_prob = zeros(L-1)
Gen_prob = zeros(G-1)
for l in 1:(L-1)
    bin_l = Line.line_bins[l]
    bin_h = Line.line_bins[l+1]
    Line_prob[l] = find_bins_pb(bin_l, bin_h,TotalL,true,false)
end
for g in 1:(G-1)
    bin_l = Gen.gen_bins[j]
    bin_h = Gen.gen_bins[j+1]
    Gen_prob[g] = find_bins_pb(bin_l, bin_h,TotalG,false,true)
end
