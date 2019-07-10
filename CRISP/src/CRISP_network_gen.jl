# network segmenting functions for CRISP-electricity
# should also check out JuliaGraphs -> collection of pkgs for common graph analysis
using SparseArrays
using LinearAlgebra
using DataFrames
using CSV
#import ps from csv files
function import_ps(filename)
    psBusData = CSV.read("$filename\\bus.csv",allowmissing=:none)
    n = length(psBusData.id);
    psBusIndex =  sparse(psBusData.id,fill(1,n),collect(1:n));
    #psBusIndex = CSV.read("$filename\\bi.csv",allowmissing=:none)
    psBranchData = CSV.read("$filename\\branch.csv",allowmissing=:none);
    psGenData = CSV.read("$filename\\gen.csv",allowmissing=:none);
    psShuntData = CSV.read("$filename\\shunt.csv",allowmissing=:none);
    if isfile("$filename\\storage.csv") psStorageData = CSV.read("$filename\\storage.csv",allowmissing=:none);
    else psStorageData = DataFrame(bus = [], E = [], Ps = [], Emax = [], Emin = [], Psmax = [], Psmin = [], status = []); end
    mpBaseMVA =  100; # CSV.read("$filename\\baseMVA.csv")[1,1];
    #if !isempty(ps.gencost) CSV.write("$filename-gen_cost.csv",ps.gencost) end
    ## Changing types in dataframes:
    psBranchData[:,:Pf] = psBranchData[:,:Pf].*1.0;
    psBranchData[:,:Qf] = psBranchData[:,:Qf].*1.0;
    psBranchData[:,:Pt] = psBranchData[:,:Pt].*1.0;
    psBranchData[:,:Qt] = psBranchData[:,:Qt].*1.0;
    psGenData[:,:Pg] = psGenData[:,:Pg].*1.0;
    psGenData[:,:Pmax] = psGenData[:,:Pmax].*1.0;
    psShuntData[:,:P] = psShuntData[:,:P].*1.0;
    ps = PSCase(mpBaseMVA, psBusData, psBranchData, psGenData, psShuntData, psStorageData, psBusIndex);
    return ps
end

# exports ps structure to several csv files
function export_ps(ps,filename)
    if !isempty(ps.bus) CSV.write("$filename\\bus.csv",ps.bus) end
    if !isempty(ps.branch) CSV.write("$filename\\branch.csv",ps.branch) end
    if !isempty(ps.gen) CSV.write("$filename\\gen.csv",ps.gen) end
    if !isempty(ps.shunt) CSV.write("$filename\\shunt.csv",ps.shunt) end
    if !isempty(ps.baseMVA) CSV.write("$filename\\baseMVA.csv",DataFrame(base_MVA = ps.baseMVA)) end
    if !isempty(ps.storage) CSV.write("$filename\\storage.csv",ps.storage) end
    #if !isempty(ps.bi)
        #n = length(ps.bus.id);
        #bi = sparse(ps.bus.id,fill(1,n),collect(1:n));
        #CSV.write("$filename\\bi.csv",bi)
    #end
end

function find_subgraphs(ps)
    # usage: [graphNos,nSubGraphs,linkNos] = find_subgraphs(nodes_A,links)
    #  where nodes is an n x 1 list of node numbers and links is
    #  m x (2+) list of edges (from, to)
    # The return value is a n x 1 vector of sub-graph numbers
    subgraphs = zeros(length(ps.bus[:id]))

    n = size(ps.bus,1) # the number of buses
    bi = sparse(ps.bus.id,fill(1,n),collect(1:n))
    stats = (ps.branch[:status].==1);
    links = [bi[ps.branch[stats,:f]] bi[ps.branch[stats,:t]]]

    n = length(ps.bus[:id]);
    m = length(ps.branch[stats,:f]);

    A = zeros(n,n);
    for i = 1:m
        A[links[i,1],links[i,2]] = 1
        A[links[i,2],links[i,1]] = 1
    end
    A = A+Matrix{Float64}(I,n,n)
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
    return subgraphs = Int64.(graphNos)
end

struct Island_ps
    bus::Array
    branch::Array
    shunt::Array
    gen::Array
    storage::Array
end

function build_islands(subgraph,ps)
    N = Int64(findmax(subgraph)[1]);
    ps_islands = Array{Island_ps}(undef,N);
    for jj = 1:N
        nodes = subgraph.==jj;
        buses = ps.bus[nodes,:id];
        gen = falses(length(ps.gen[:bus]));
        shunt = falses(length(ps.shunt[:bus]));
        storage = falses(length(ps.storage[:bus]));
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
        for st = 1:length(ps.storage[:bus])
            if sum(ps.storage[st,:bus].==buses)!=0
                storage[st] = true;
            end
        end
        for l = 1:length(ps.branch[:f])
            if ps.branch[l,:status]!=0
                if sum(ps.branch[l,:f].==buses)!=0 && sum(ps.branch[l,:t].==buses)!=0
                  branch[l]=true;
                end
            end
        end
        ps_islands[jj] = Island_ps(nodes,branch,shunt,gen,storage);
    end
    return ps_islands
end

mutable struct PSCase
    baseMVA::Int64
    bus::DataFrame
    branch::DataFrame
    gen::DataFrame
    shunt::DataFrame
    storage::DataFrame
    bi::SparseMatrixCSC{Int64,Int64} # TODO: figure out how to make this the same as Pavan's PSCase
end

function ps_subset(ps,ps_island)
    mpBaseMVA = ps.baseMVA;
    psBusData = ps.bus[ps_island.bus,:];
    psBranchData = ps.branch[ps_island.branch,:];
    psGenData = ps.gen[ps_island.gen,:];
    psShuntData = ps.shunt[ps_island.shunt,:];
    psStorageData = ps.storage[ps_island.storage,:];
    n = length(psBusData.id);
    psBusIndex =  sparse(psBusData.id,fill(1,n),collect(1:n));#DataFrame(bi=ps.bi.bi[ps_island.bus]);
    psi = PSCase(mpBaseMVA, psBusData, psBranchData, psGenData, psShuntData, psStorageData, psBusIndex);
    return psi
end

function add_changes!(ps,psi,ps_island);
    ps.gen[ps_island.gen,:Pg] = psi.gen.Pg
    ps.shunt[ps_island.shunt,:P] = psi.shunt.P
    ps.storage[ps_island.storage,:P] = psi.storage.P
    ps.storage[ps_island.storage,:E] = psi.storage.E
end

function ps_visualize(ps,filename)
    data = []
    nbu = size(ps.bus,1)
    nbr = size(ps.branch,1)
    # initialize bus sizes
    bus_sizes = zeros(nbu,1)

    for i = 1:nbr
        # TODO: Convert this using the bus index
        # find the end points
        f = ps.branch.f[i]
        t = ps.branch.t[i]
        fromX = ps.bus.locX[f]
        toX = ps.bus.locX[t]
        fromY = ps.bus.locY[f]
        toY = ps.bus.locY[t]
        # choose the line thickness
        br_size = sqrt(abs(ps.branch.Pf[i])) * 4;
        # update the bus size
        bus_sizes[f] = max(br_size,bus_sizes[f])
        bus_sizes[t] = max(br_size,bus_sizes[t])
        # choose the line color
        # TODO: compute the line color here
        # plot
        line = scatter(;x=[fromX,toX],y=[fromY,toY],mode="lines",line_width=br_size)
        if i==1
            data = [line]
        else
            push!(data,line)
        end
    end
    # plot the buses
    nbu = size(ps.bus,1)
    for i=1:nbu
        # decide on the marker size
        s = bus_sizes[i]
        bus = scatter(;x=[ps.bus.locX[i]],y=[ps.bus.locY[i]],mode="markers",marker_size=s)
        push!(data,bus)
    end

    println(data)
    out = plot(data,Layout(title="A power system"))
    f = "$filename.pdf"
    savefig(out,f)
    return f
end