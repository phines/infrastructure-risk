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
    A = A+eye(length(ps.bus[:id]))

    grNo = 1;
    graphNos = zeros(n);
    linkNos_int = zeros(m_int);
    next = 1;
    while ~isempty(next)
        included = false(n,1);
        included(next) = true;
        oldLen = 0;
        while sum(included) ~= oldLen
            oldLen = sum(included);
            [Ai,~] = find(A(:,included));
            included(Ai) = true;
      end
      graphNos(included) = grNo;
      grNo = grNo+1;
      next = find(graphNos==0,1);
end
    return subgraphs
end

function subsetps(subgraphs,ps)

    return ps_struct
end
