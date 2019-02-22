# network segmenting functions for CRISP-electricity
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

function subsetps(ps,subgraph)
    
    return ps_subgraph
end
