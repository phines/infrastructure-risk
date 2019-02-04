# network segmenting functions for CRISP-electricity
function find_subgraphs(ps)
    # usage: [graphNos,nSubGraphs,linkNos] = find_subgraphs(nodes_A,links)
    #  where nodes is an n x 1 list of node numbers and links is
    #  m x (2+) list of edges (from, to)
    # The return value is a n x 1 vector of sub-graph numbers
    subgraphs = zeros(length(ps.bus[:id]))
    stats = (ps.branch[:status].==1);
    links = [ps.branch[stats,:f] ps.branch[stats,:t]]
    A = zeros(length(ps.bus[:id]),length(ps.bus[:id]));
    A[links[1,:],links[2,:]] = 1

    return subgraphs
end

function subsetps(subgraphs,ps)

    return ps_struct
end
