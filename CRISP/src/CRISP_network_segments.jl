# network segmenting functions for CRISP-electricity
# should also check out JuliaGraphs -> collection of pkgs for common graph analysis
using SparseArrays
using LinearAlgebra
using DataFrames

function find_subgraphs(ps)
    # usage: [graphNos,nSubGraphs,linkNos] = find_subgraphs(nodes_A,links)
    #  where nodes is an n x 1 list of node numbers and links is
    #  m x (2+) list of edges (from, to)
    # The return value is a n x 1 vector of sub-graph numbers
    subgraphs = zeros(length(ps.bus[:id]))

    stats = (ps.branch[:status].==1);
    links = [ps.branch[stats,:f] ps.branch[stats,:t]]

    n = length(ps.bus[:id]);
    m = length(links[:,1]);

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
        included = zeros(n).==1;
        included[index] = true;
        oldLen = 0;
        while sum(included) ~= oldLen
            oldLen = sum(included);
            Ai = A[index,:].==1;
            included[Ai] .= true;
        end
        graphNos[included] .= grNo;
        grNo = grNo+1;
        all_indices = (LinearIndices(graphNos))[graphNos.==0]
        index = all_indices[1];
    end
    return subgraphs = graphNos
end

<<<<<<< HEAD
function subsetps(ps,subgraph)
    
    return ps_subgraph
=======
function build_islands(subgraph,ps)
    N = max(subgraphs);
    struct island
        bus::Array
        branch::Array
        shunt::Array
        gen::Array
    end
    ps_ilands = Array{iland}(N);
    for jj = 1:N
        nodes = subgraph.==jj;
        gen = falses(length(ps.gen[:bus]));
        shunt = falses(length(ps.shunt[:bus]));
        branch = falses(length(ps.branch[:f]));
        for g = 1:length(ps.gen[:bus])
            if sum(ps.gen[:bus].==nodes)
                gen[g] = true;
            end
        end
        for s = 1:length(ps.shunt[:bus])
            if sum(ps.shunt[:bus].==nodes)
                shunt[s] = true;
            end
        end
        for l = 1:length(ps.branch[:f])
            if ps.branch[:status]
                if ps.branch[f].==l
                  branch[l]=true;
                end
                if ps.branch[t].==l
                  branch[l]=true;
                end
            end
        end
        ps_islands[jj] = island(nodes,branch,shunt,gen)
    end
    return ps_islands
>>>>>>> 8a61a1e4eca49e36cac2202c996445b8ad5870de
end
