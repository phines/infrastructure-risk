# network segmenting functions for CRISP-electricity
# should also check out JuliaGraphs -> collection of pkgs for common graph analysis
using SparseArrays
using LinearAlgebra
using DataFrames
using CSV

function reduce_network!(ps,deg)
    tolerance = 10^(-6)
    nl = length(ps.shunt.bus)
    mismatch = sum(ps.shunt.P)-sum(ps.gen.Pg)-sum(ps.hvdc.P)
    Index = 1:length(ps.bus.id)
    m = Index[deg .== 1]
    mis_set_P = -mismatch/(length(m))
    pf = 0.95 #assuming powerfactor of 0.85
    mis_set_Q = sin(acos(pf))*(mis_set_P/pf) #assuming pf of 0.95 for loads
    j=0;
    ind_l = Int64.(zeros(length(m)))
    for i in m
        Index_lines = 1:length(ps.branch.f)
        Index_bus = 1:length(ps.bus.id)
        j=j+1;
        Fr = (ps.branch.f .== ps.bus.id[i])
        To = (ps.branch.t .== ps.bus.id[i])
        if sum(Fr)==0
            b = ps.branch.f[To][1]
            name = ps.bus.name[ps.bus.id .== b][1]
            ps.shunt.bus[ps.shunt.bus .== ps.bus.id[i]] .= b
            ps.shunt.name[ps.shunt.bus .== ps.bus.id[i]] .= name
            ps.gen.bus[ps.gen.bus .== ps.bus.id[i]] .= b
            ps.gen.name[ps.gen.bus .== ps.bus.id[i]] .= name
            if mismatch < -tolerance
                l = DataFrame(id = nl+j,bus=b, name =name*"_mismatch", status = 1.0, P = mis_set_P, Q=mis_set_Q);
                append!(ps.shunt,l)
            end
            ind_l[j] = Index_lines[ps.branch.t .== ps.bus.id[i]][1];
        else
            b = ps.branch.t[Fr][1]
            name = ps.bus.name[ps.bus.id .== b][1]
            ps.shunt.bus[ps.shunt.bus .== ps.bus.id[i]] .= b
            ps.shunt.name[ps.shunt.bus .== ps.bus.id[i]] .= name
            ps.gen.bus[ps.gen.bus .== ps.bus.id[i]] .= b
            ps.gen.name[ps.gen.bus .== ps.bus.id[i]] .= name
            if mismatch < -tolerance
                l = DataFrame(id = nl+j,bus=b, name =name*"_mismatch", status = 1.0, P = mis_set_P, Q=mis_set_Q);
                append!(ps.shunt,l)
            end
            ind_l[j] = Index_lines[ps.branch.f .== ps.bus.id[i]][1];
        end
    end
    sort!(ind_l)
    deleterows!(ps.branch,ind_l)
    deleterows!(ps.bus,m)
    return ps
end

function find_node_degree(nodes,f,t)
    n = length(nodes)
    degree = zeros(n)
    for n in 1:n
        degree[n] = sum(nodes[n] .== f) + sum(nodes[n] .== t);
    end
    return degree
end

function choose_gens_black_start!(ps,fraction, sizethreshold)
    ng = length(ps.gen.bus)
    ps.gen[!,:black_start] = falses(ng)
    st = ps.gen.Pmax .<= sizethreshold
    fr = Int64(round(sum(st)*fraction))
    if sum(st) >= 1
        ps.gen.black_start[st][rand(rng,(1:sum(st)),fr)] .= true
    end
    ps.gen.service_load[ps.gen.black_start .== true] .= 0
    return ps
end

# sets generator states for modeling black start
function set_gen_states!(ps)
    for g in 1:length(ps.gen.Pg)
        if ps.gen.Pg[g] > 0
            ps.gen.state[g] = On;
        elseif ps.gen.status[g] == 0
            ps.gen.state[g] = OutOfOpperation;
            ps.gen.service_load[g] = 0;
        else
            ps.gen.state[g] = Off
            ps.gen.service_load[g] = 0;
        end
    end
    return ps
end

mutable struct PSCase
    baseMVA::Int64
    bus::DataFrame
    branch::DataFrame
    gen::DataFrame
    shunt::DataFrame
    storage::DataFrame
    hvdc::DataFrame
    bi::SparseMatrixCSC{Int64,Int64}
end

#import ps from csv files
function import_ps(filename)
    psBusData = CSV.File("$filename/bus.csv")  |> DataFrame
    n = length(psBusData.id);
    psBusIndex =  sparse(psBusData.id,fill(1,n),collect(1:n));
    #psBusIndex = CSV.read("$filename/bi.csv",allowmissing=:none)
    psBranchData = CSV.File("$filename/branch.csv")  |> DataFrame;
    if occursin("bs",filename)
        psGenData = CSV.File("$filename/gen.csv",
               types = [Int64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64,
               Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64,
               Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, String,
               Float64, Float64, Float64, Float64, String, Float64, String, String, String, String, Float64,
               Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64,
               String, Float64, String, String, String, String,Float64,Float64,Bool]) |> DataFrame;
    elseif occursin("case73_",filename)
        psGenData = CSV.File("$filename/gen.csv",
               types = [Int64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64,
               Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64,
               Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, String,
               Float64, Float64, Float64, Float64, String, Float64, String, String, String, String, Float64,
               Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64,
               String, Float64, String, String, String, String]) |> DataFrame;
    else
        psGenData = CSV.File("$filename/gen.csv") |> DataFrame;
    end
    psShuntData = CSV.File("$filename/shunt.csv")  |> DataFrame;
    mpBaseMVA =  100; # CSV.read("$filename/baseMVA.csv")[1,1];
    #if !isempty(ps.gencost) CSV.write("$filename-gen_cost.csv",ps.gencost) end
    ## Changing types in dataframes:
    if occursin("NE_NY_",filename)
        psBranchData[!,:Pf] = zeros(length(psBranchData.f)) .* 1.0;
        psBranchData[!,:Qf] = zeros(length(psBranchData.f)) .* 1.0;
        psBranchData[!,:Pt] = zeros(length(psBranchData.f)) .* 1.0;
        psBranchData[!,:Qt] = zeros(length(psBranchData.f)) .* 1.0;
    else
        psBranchData.Pf = psBranchData.Pf .* 1.0;
        psBranchData.Qf = psBranchData.Qf .* 1.0;
        psBranchData.Pt = psBranchData.Pt .* 1.0;
        psBranchData.Qt = psBranchData.Qt .* 1.0;
        psGenData.Pg = psGenData.Pg .* 1.0;
        psGenData.Qg = psGenData.Qg .* 1.0;
    end
    psGenData.Pmax = psGenData.Pmax .* 1.0;
    psGenData.Pmin = psGenData.Pmin .* 1.0;
    psShuntData.P = psShuntData.P .* 1.0;
    psShuntData.status = psShuntData.status .* 1.0;
    if isfile("$filename/storage.csv") psStorageData = CSV.File("$filename/storage.csv")   |> DataFrame;
        psStorageData.Ps = psStorageData.Ps .* 1.0;
        psStorageData.E = psStorageData.E .* 1.0;
        psStorageData.Psmax = psStorageData.Psmax .* 1.0;
        psStorageData.Psmin = psStorageData.Psmin .* 1.0;
        psStorageData.Emax = psStorageData.Emax .* 1.0;
        psStorageData.Emin = psStorageData.Emin .* 1.0;
    else psStorageData = DataFrame(bus = Int64[], E = Float64[], Ps = Float64[], Emax = Float64[], Emin = Float64[], Psmax = Float64[], Psmin = Float64[], Efficiency = Float64[], status = Int64[]);
    end
    if isfile("$filename/hvdc.csv") psHVDCData = CSV.File("$filename/hvdc.csv")   |> DataFrame;
    else psHVDCData = DataFrame(id = Int64[], bus = Int64[], name = String[], status = Int64[], P = Float64[], Q = Float64[], Pmax = Float64[], Qmax = Float64[]);
    end
    # add status column if there is not one
    if (sum(:status .== names(psShuntData)) < 1)
        psShuntData[!,:status] = ones(size(psShuntData,1))
    end
    if (sum(:status .== names(psGenData)) < 1)
        psGenData[!,:status] = ones(size(psGenData,1))
    end
    if (sum(:status .== names(psBranchData)) < 1)
        psBranchData[!,:status] = ones(size(psBranchData,1))
    end
    if (sum(:status .== names(psStorageData)) < 1)
        psStorageData[!,:status] = ones(size(psStorageData,1))
    end
    ps = PSCase(mpBaseMVA, psBusData, psBranchData, psGenData, psShuntData, psStorageData, psHVDCData, psBusIndex);
    return ps
end

# exports ps structure to several csv files
function export_ps(ps,filename)
    if !isempty(ps.bus) CSV.write("$filename/bus.csv",ps.bus) end
    if !isempty(ps.branch) CSV.write("$filename/branch.csv",ps.branch) end
    if !isempty(ps.gen) CSV.write("$filename/gen.csv",ps.gen) end
    if !isempty(ps.shunt) CSV.write("$filename/shunt.csv",ps.shunt) end
    if !isempty(ps.baseMVA) CSV.write("$filename/baseMVA.csv",DataFrame(base_MVA = ps.baseMVA)) end
    if !isempty(ps.storage) CSV.write("$filename/storage.csv",ps.storage) end
    if !isempty(ps.hvdc) CSV.write("$filename/hvdc.csv",ps.hvdc) end
end

function find_subgraphs(ps)
    # usage: [graphNos,nSubGraphs,linkNos] = find_subgraphs(nodes_A,links)
    #  where nodes is an n x 1 list of node numbers and links is
    #  m x (2+) list of edges (from, to)
    # The return value is a n x 1 vector of sub-graph numbers
    subgraphs = zeros(length(ps.bus.id))

    n = size(ps.bus,1) # the number of buses
    bi = sparse(ps.bus.id,fill(1,n),collect(1:n))
    stats = (ps.branch.status .== 1);
    links = [bi[ps.branch.f[stats]] bi[ps.branch.t[stats]]]

    n = length(ps.bus.id);
    m = length(ps.branch.f[stats]);

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
    hvdc::Array
end

function build_islands(subgraph,ps)
    N = Int64(findmax(subgraph)[1]);
    ps_islands = Array{Island_ps}(undef,N);
    for jj = 1:N
        nodes = subgraph.==jj;
        buses = ps.bus[nodes,:id];
        gen = falses(length(ps.gen.bus));
        shunt = falses(length(ps.shunt.bus));
        storage = falses(length(ps.storage.bus));
        branch = falses(length(ps.branch.f));
        hvdc = falses(length(ps.hvdc.bus))
        for g = 1:length(ps.gen.bus)
            if sum(ps.gen.bus[g] .== buses)!=0
                gen[g] = true;
            end
        end
        for h in 1:length(ps.hvdc.bus)
            if sum(ps.hvdc.bus[h] .== buses)!=0
                hvdc[h] = true
            end
        end
        for s = 1:length(ps.shunt.bus)
            if sum(ps.shunt.bus[s].==buses)!=0
                shunt[s] = true;
            end
        end
        for st = 1:length(ps.storage.bus)
            if sum(ps.storage.bus[st] .== buses)!=0
                storage[st] = true;
            end
        end
        for l = 1:length(ps.branch.f)
            if ps.branch.status[l] != 0
                if (sum(ps.branch.f[l] .== buses) != 0) && (sum(ps.branch.t[l] .== buses) != 0)
                  branch[l]=true;
                end
            end
        end
        ps_islands[jj] = Island_ps(nodes,branch,shunt,gen,storage,hvdc);
    end
    return ps_islands
end

function find_diameter(ps)
    nodes = ps.bus.id
    n = length(nodes)
    edges = [ps.branch.f ps.branch.t]
    diameter = zeros(n)
    subset = Array{Int64,1}(undef,1)
    for s in 1:n
        k = 1
        subset[1] = nodes[s]
        n1 = [0]
        n2 = [1 1]
        while length(n1) < length(n2)
            diameter[s] = k
            n1 = find_neighbors(nodes, edges, subset; K=k)
            n2 = find_neighbors(nodes, edges, subset; K=(k+1))
            k = k+1
        end
    end
    diameter = maximum(diameter)
    return Int64(diameter)
end

function find_neighbors(nodes, edges, subset; K=1)
    #finds all neighbors of degree <=K for subset node
    n = length(nodes)
    neighbors = falses(n)
    F = edges[:,1]
    T = edges[:,2]
    ei = sparse(nodes, repeat([1],inner=n), 1:n, maximum(nodes), 1)
    k  = intersect(findall(.!isnothing.(indexin(edges[:,1],nodes))), findall(.!isnothing.(indexin(edges[:,2],nodes)))) # link indices for which (from_node,to_node) belong to nodes
	f  = ei[edges[k,1]]
	t  = ei[edges[k,2]]
    A = sparse([f;t],[t;f],1,n,n)
    (A[A .>= 1] .= 1)
	(A += SparseMatrixCSC{Float64}(I,n,n))
    neighbors[ei[subset]] .= true
	for k in 1:K
		(i, _, _) = findnz(A[:,neighbors])
		neighbors[i] .= true
	end
	neighbors[ei[subset]] .= false
    return nodes[neighbors]
end

function index_in(a, b)
	i = indexin(a, b)
	j = .!isnothing.(i)

	return i[j], j
end

function find_lines_n_hops(ps,lines_status,hop)
    subset = Array{Int64,1}(undef,2)
    nodes = ps.bus.id
    edges = [ps.branch.f ps.branch.t]
    node_neighbor = [];
    new_line = [];
    counter = 0
    while (isempty(node_neighbor) || isempty(new_line)) && (counter <= (sum(lines_status .== 0)*3))
        counter = counter + 1
        i = rand(rng,1:sum(lines_status .== 0))
        subset = edges[lines_status .== 0,:][i,:]
        if hop > 1
            node_neighbor =  find_neighbors(nodes, edges, subset; K=hop)
            node_neighbor_less =  find_neighbors(nodes, edges, subset; K=(hop-1))
            new_line = find_line_n_neighbors(edges,node_neighbor,node_neighbor_less,subset,lines_status)
        else
            node_neighbor =  find_neighbors(nodes, edges, subset; K=hop)
            new_line = find_line_neighbors(edges,node_neighbor,subset,lines_status)
        end
    end
    return new_line
end

function find_line_neighbors(edges,node_neighbor,subset,lines_status)
    lst = falses(length(lines_status))
    lst[lines_status .== 1] .= true
    nl = length(lines_status)
    n = length(node_neighbor)
    Index = 1:nl
    cascade = Index[lst]
    possible_lines = falses(nl)
    for h in 1:n
        find_edges_t = edges[:,1] .== node_neighbor[h]
        possible_lines[find_edges_t] .= true
        find_edges_f = edges[:,2] .== node_neighbor[h]
        possible_lines[find_edges_f] .= true
    end
    #make sure there's at least one possible line
    if sum(possible_lines .== 1) == 0
        return new_line = []
    end
    #remove lines that are parallel to subset
    n11 = edges[:,1] .== subset[1]
    possible_lines[n11] .= false
    n21 = edges[:,2] .== subset[1]
    possible_lines[n21] .= false
    n12 = edges[:,1] .== subset[2]
    possible_lines[n12] .= false
    n22 = edges[:,2] .== subset[2]
    possible_lines[n22] .= false
    #remove lines that are already in cascade
    possible_lines[lines_status .== 0] .= false
    if sum(possible_lines .== 1) == 0
        return new_lines = []
    else
        #pick line
        line_nums = Index[possible_lines]
        p = rand(rng,1:length(line_nums))
        new_line = line_nums[p]
    end
    return new_line
end

function find_line_n_neighbors(edges,node_neighbor,node_neighbor_less,subset,lines_status)
    lst = falses(length(lines_status))
    lst[lines_status .== 1] .= true
    neighbors = []
    less = subset[1]
    less = [less subset[2]]
    for n in 1:length(node_neighbor)
        if sum(node_neighbor[n] .== node_neighbor_less) > 0
            less = [less node_neighbor[n]]
        else
            if isempty(neighbors)
                neighbors = node_neighbor[n]
            else
                neighbors = [neighbors node_neighbor[n]]
            end
        end
    end
    nl = length(lines_status)
    n = length(node_neighbor)
    Index = 1:nl
    cascade = Index[lst]
    possible_lines = falses(nl)
    for h in 1:n
        find_edges_t = edges[:,1] .== node_neighbor[h]
        possible_lines[find_edges_t] .= true
        find_edges_f = edges[:,2] .== node_neighbor[h]
        possible_lines[find_edges_f] .= true
    end
    #make sure there's at least one possible line
    if sum(possible_lines .== 1) == 0
        return new_line = []
    end
    #remove lines that are already in cascade
    possible_lines[lines_status .== 0] .= false
    if sum(possible_lines .== 1) == 0
        return new_lines = []
    end
    #remove closer neighbors including in the subset
    for m in 1:length(less)
        n11 = edges[:,1] .== less[m]
        possible_lines[n11] .= false
        n21 = edges[:,2] .== less[m]
        possible_lines[n21] .= false
    end
    if sum(possible_lines .== 1) == 0
        return new_lines = []
    else
        #pick line
        line_nums = Index[possible_lines]
        p = rand(rng,1:length(line_nums))
        new_line = line_nums[p]
    end
    return new_line
end

function ps_subset(ps,ps_island)
    mpBaseMVA = ps.baseMVA;
    psBusData = ps.bus[ps_island.bus,:];
    psBranchData = ps.branch[ps_island.branch,:];
    psGenData = ps.gen[ps_island.gen,:];
    psShuntData = ps.shunt[ps_island.shunt,:];
    psStorageData = ps.storage[ps_island.storage,:];
    psHVDCData = ps.hvdc[ps_island.hvdc,:];
    n = length(psBusData.id);
    psBusIndex =  sparse(psBusData.id,fill(1,n),collect(1:n));#DataFrame(bi=ps.bi.bi[ps_island.bus]);
    psi = PSCase(mpBaseMVA, psBusData, psBranchData, psGenData, psShuntData, psStorageData, psHVDCData, psBusIndex);
    return psi
end

function add_changes!(ps,psi,ps_island);
    ps.gen.Pg[ps_island.gen] = psi.gen.Pg
    if sum(names(ps.gen).==:time_on) == 0
        ps.gen.time_off[ps_island.gen] = psi.gen.time_off
        ps.gen.time_on[ps_island.gen] = psi.gen.time_on
    end
    ps.shunt.status[ps_island.shunt] = psi.shunt.status
    ps.storage.Ps[ps_island.storage] = psi.storage.Ps
    ps.storage.E[ps_island.storage] = psi.storage.E
end

#import ps from csv files
function import_ps0(filename)
    #psBusIndex = CSV.File("$filename/bi.csv") |> DataFrame
    psBusData = CSV.File("$filename/bus.csv") |> DataFrame
    psBranchData = CSV.File("$filename/branch.csv") |> DataFrame
    psGenData = CSV.File("$filename/gen.csv") |> DataFrame
    psShuntData = CSV.File("$filename/shunt.csv") |> DataFrame
    mpBaseMVA =  100; # CSV.read("$filename/baseMVA.csv")[1,1];
    #if !isempty(ps.gencost) CSV.write("$filename-gen_cost.csv",ps.gencost) end
    ## Changing types in dataframes:
    psBranchData.Pf = psBranchData.Pf.*1.0;
    psBranchData.Qf = psBranchData.Qf.*1.0;
    psBranchData.Pt = psBranchData.Pt.*1.0;
    psBranchData.Qt = psBranchData.Qt.*1.0;
    psGenData.Pg = psGenData.Pg.*1.0;
    psGenData.Pmax = psGenData.Pmax.*1.0;
    psShuntData.P = psShuntData.P.*1.0;
    n = length(psBusData.id);
    bi = sparse(psBusData.id,fill(1,n),collect(1:n));
    psBusIndex = bi;
    ps = PSCase0(mpBaseMVA, psBusData, psBranchData, psGenData, psShuntData, psBusIndex);
    return ps
end

# exports ps structure to several csv files
function export_ps0(ps,filename)
    if !isempty(ps.bus) CSV.write("$filename/bus.csv",ps.bus) end
    if !isempty(ps.branch) CSV.write("$filename/branch.csv",ps.branch) end
    if !isempty(ps.gen) CSV.write("$filename/gen.csv",ps.gen) end
    if !isempty(ps.shunt) CSV.write("$filename/shunt.csv",ps.shunt) end
    if !isempty(ps.baseMVA) CSV.write("$filename/baseMVA.csv",DataFrame(base_MVA = ps.baseMVA)) end
    if !isempty(ps.bi)
        #n = length(ps.bus.id);
        #bi = sparse(ps.bus.id,fill(1,n),collect(1:n));
        CSV.write("$filename/bi.csv",ps.bi)
    end
end

function find_subgraphs0(ps)
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

struct Island_ps0
    bus::Array
    branch::Array
    shunt::Array
    gen::Array
end

function build_islands0(subgraph,ps)
    N = Int64(findmax(subgraph)[1]);
    ps_islands = Array{Island_ps}(undef,N);
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
        ps_islands[jj] = Island_ps0(nodes,branch,shunt,gen);
    end
    return ps_islands
end

mutable struct PSCase0
    baseMVA::Int64
    bus::DataFrame
    branch::DataFrame
    gen::DataFrame
    shunt::DataFrame
    bi::DataFrame # not quite right, Pavan uses SparseMatrixCSC{Int64,Int64} # TODO: figure out how to make this the same as Pavan's PSCase
end

function ps_subset0(ps,ps_island)
    mpBaseMVA = ps.baseMVA;
    psBusData = ps.bus[ps_island.bus,:];
    psBranchData = ps.branch[ps_island.branch,:];
    psGenData = ps.gen[ps_island.gen,:];
    psShuntData = ps.shunt[ps_island.shunt,:];
    n = length(psBusData.id);
    psBusIndex =  sparse(psBusData.id,fill(1,n),collect(1:n));#psBusIndex = DataFrame(bi=ps.bi.bi[ps_island.bus]);
    psi = PSCase0(mpBaseMVA, psBusData, psBranchData, psGenData, psShuntData, psBusIndex);
    return psi
end

function add_changes0!(ps,psi,ps_island);
    ps.gen[ps_island.gen,:Pg] = psi.gen.Pg
    ps.shunt[ps_island.shunt,:P] = psi.shunt.P
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
