# exports ps structure to several csv files
using CSV; using DataFrames
function export_ps(ps,filename)
    if !isempty(ps.bus) CSV.write("$filename-bus.csv",ps.bus) end
    if !isempty(ps.branch) CSV.write("$filename-branch.csv",ps.branch) end
    if !isempty(ps.gen) CSV.write("$filename-gen.csv",ps.gen) end
    if !isempty(ps.shunt) CSV.write("$filename-shunt.csv",ps.shunt) end
    if !isempty(ps.baseMVA) CSV.write("$filename-baseMVA.csv",DataFrame(base_MVA = ps.baseMVA)) end
    if !isempty(ps.gencost) CSV.write("$filename-gen_cost.csv",ps.gencost) end
end
