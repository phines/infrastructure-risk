# network segmenting functions for CRISP-electricity
# should also check out JuliaGraphs -> collection of pkgs for common graph analysis
using SparseArrays
using LinearAlgebra
using DataFrames
using CSV
#import ps from csv files
function import_ps(filename)
    psBusIndex = CSV.read("$filename\\bus.csv")
    psBusData = CSV.read("$filename\\bus.csv")
    psBranchData = CSV.read("$filename\\branch.csv");
    psGenData = CSV.read("$filename\\gen.csv");
    psShuntData = CSV.read CSV.write("$filename\\shunt.csv");
    mpBaseMVA = CSV.write("$filename\\baseMVA.csv")[1,1];
    #if !isempty(ps.gencost) CSV.write("$filename-gen_cost.csv",ps.gencost) end
    ps = PSCase(mpBaseMVA, psBusData, psBranchData, psGenData, psGenCostData, psShuntData);
    return ps
end

# exports ps structure to several csv files
function export_ps(ps,filename)
    if !isempty(ps.bus) CSV.write("$filename-bus.csv",ps.bus) end
    if !isempty(ps.branch) CSV.write("$filename-branch.csv",ps.branch) end
    if !isempty(ps.gen) CSV.write("$filename-gen.csv",ps.gen) end
    if !isempty(ps.shunt) CSV.write("$filename-shunt.csv",ps.shunt) end
    if !isempty(ps.baseMVA) CSV.write("$filename-baseMVA.csv",DataFrame(base_MVA = ps.baseMVA)) end
    if !isempty(ps.bi) CSV.write("$filename-bi.csv",ps.bi) end
end

function find_subgraphs(ps)
    # usage: [graphNos,nSubGraphs,linkNos] = find_subgraphs(nodes_A,links)
    #  where nodes is an n x 1 list of node numbers and links is
    #  m x (2+) list of edges (from, to)
    # The return value is a n x 1 vector of sub-graph numbers
    subgraphs = zeros(length(ps.bus[:id]))

    stats = (ps.branch[:status].==1);
    links = [ps.branch[stats,:f] ps.branch[stats,:t]]

    n = length(ps.bus[:id]);
    m = length(ps.branch[stats,:f]);

    A = zeros(n,n);
    for i = 1:m
        A[links[i,1],links[i,2]] = 1
        A[links[i,2],links[i,1]] = 1
    end
    A = A+eye(n)
    m_int = m; #setting the internal links to all the links
    grNo = 1;
    graphNos = zeros(n);
    linkNos_int = zeros(m_int);
    index = 1;
    while ~isempty(index)
        included = falses(n);
        included[index] = true;
        oldLen = 0;
        while sum(included) != oldLen
            oldLen = sum(included);
            for jj = 1:n
                ii = included[jj];
                if ii
                    Ai = A[:,jj].==1;
                    included[Ai] .= true;
                end
            end
        end
        graphNos[included] .= grNo;
        grNo = grNo+1;
        all_indices = (LinearIndices(graphNos))[graphNos.==0];
        if !isempty(all_indices)
            index = all_indices[1];
        else
            index = (LinearIndices(graphNos))[graphNos.==0];
        end
    end
    return subgraphs = graphNos
end

struct Island_ps
    bus::Array
    branch::Array
    shunt::Array
    gen::Array
end

function build_islands(subgraph,ps)
    N = Int64(findmax(subgraph)[1]);
    ps_islands = Array{Island_ps}(N);
    for jj = 1:N
        nodes = subgraph.==jj;
        buses = ps.bus[nodes,:id];
        gen = falses(length(ps.gen[:bus]));
        shunt = falses(length(ps.shunt[:bus]));
        branch = falses(length(ps.branch[:f]));
        for g = 1:length(ps.gen[:bus])
            if sum(ps.gen[g,:bus].==buses)!=0
                gen[g] = true;
            end
        end
        for s = 1:length(ps.shunt[:bus])
            if sum(ps.shunt[s,:bus].==buses)!=0
                shunt[s] = true;
            end
        end
        for l = 1:length(ps.branch[:f])
            if ps.branch[l,:status]!=0
                if sum(ps.branch[l,:f].==buses)!=0 && sum(ps.branch[l,:t].==buses)!=0
                  branch[l]=true;
                end
            end
        end
        ps_islands[jj] = Island_ps(nodes,branch,shunt,gen)
    end
    return ps_islands
end

mutable struct PSCase
    baseMVA::Int64
    bus::DataFrame
    branch::DataFrame
    gen::DataFrame
    shunt::DataFrame
    bi::DataFrame
end

function ps_subset(ps,ps_island)
    mpBaseMVA = ps.baseMVA;
    psBusData = ps.bus[ps_island.bus,:];
    psBranchData = ps.branch[ps_island.branch,:];
    psGenData = ps.gen[ps_island.gen,:];
    psShuntData = ps.shunt[ps_island.shunt,:];
    psBusIndex = ps.bi[ps_island.bi,:];
    psi = PSCase(mpBaseMVA, psBusData, psBranchData, psGenData, psBusIndex, psShuntData);
    return psi
end
