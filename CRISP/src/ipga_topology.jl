#=
The MIT License (MIT)
Copyright 2019, Pavan Racherla (pracherla9@gmail.com) and contributors.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
=#


__precompile__()

module ipga_topology

using SparseArrays, LinearAlgebra
include("ipga_miscellaneous.jl")


@doc """`(A, ei) = adjacency_matrix(nodes, links; binary=true, selfie=false)`

Construct the adjacency matrix (A) for a graph.

Input-\n
nodes:  1-column array\n
links:  2-column array [fromNodes toNodes]\n
binary: if true, then A(i,j)=0/1\n
selfie: if true, then A(i,i)=1

Output-\n
A:  adjacency matrix\n
ei: external=>internal nodes map""" ->
function adjacency_matrix(nodes::Any, links::Array{Int64,2}; binary::Bool=true, selfie::Bool=false)
	ln = length(nodes)
	ei = sparse(nodes, repeat([1],inner=ln), 1:ln, maximum(nodes), 1)
	k  = intersect(find_in(links[:,1],nodes), find_in(links[:,2],nodes)) # link indices for which (from_node,to_node) belong to nodes
	f  = ei[links[k,1]]
	t  = ei[links[k,2]]

	A = sparse([f;t], [t;f], 1, ln, ln)

	if binary
		A[A .>= 1] .= 1
	end

	if selfie
		A += sp_eye(ln)
	end

	return A, ei
end


@doc """`(subgraphs, lSubgraph) = find_subgraphs(nodes, links)`

Find the subgraphs, if any, of a graph.

Input-\n
nodes: 1-column array\n
links: 2-column array [fromNodes toNodes]

Output-\n
subgraphs: subgraph IDs of nodes\n
lSubgraph: bit-vector selecting the largest subgraph""" ->
function find_subgraphs(nodes::Any, links::Array{Int64,2})
	N = length(nodes)

	(A, _)  = adjacency_matrix(nodes, links; selfie=true)

	subgraphs = zeros(Int64, N) # subgraph ids (N*1 Int64 array)
	n = 1
	m = 1
	while m != nothing
		j = repeat([false], inner=N)
		j[m] = true
		summa = 0
		while sum(j) != summa
			summa = sum(j)
			(i, _, _)  = findnz(A[:,j])
			j[i] .= true
		end
		subgraphs[j] .= n
		n = n+1
		m = findfirst(subgraphs .== 0)
	end

	Nsg = maximum(subgraphs)
	nsg = ones(Int64, Nsg)
	for n=1:Nsg
		nsg[n] = sum(n .== subgraphs)
	end
	(_, lsg) = findmax(nsg)
	lSubgraph = subgraphs .== lsg # the largest subgraph (N*1 Bool array)

	return subgraphs, lSubgraph
end


@doc """`neighbors = find_neighbors(nodes, links, subset; K)`

For the graph specified by `nodes` and `links`, find `neighbors` within K degrees of the given `subset` of nodes.

Input-\n
nodes:  1-column array\n
links:  2-column array [fromNodes toNodes]\n
subset, for which to find neighbors;\n
K:      degrees of separation

Output-\n
neighbors.""" ->
function find_neighbors(nodes::Any, links::Array{Int64,2}, subset::Any; K::Int64=1)
	size(subset) == () && error("subset is a scalar, but it needs to be an array!")

	(A, ei) = adjacency_matrix(nodes, links)

	neighbors = falses(length(nodes))
	neighbors[ei[subset]] .= true
	for k in 1:K
		(i, _, _) = findnz(A[:,neighbors])
		neighbors[i] .= true
	end
	neighbors[ei[subset]] .= false

	return nodes[neighbors]
end


# Export all symbols except those beginning with an underscore or in _symbols_excluded:
# (inspired by http://github.com/JuliaOpt/JuMP.jl/src/JuMP.jl)
const _symbols_excluded = [Symbol(@__MODULE__), :eval, :include]
for name in names(@__MODULE__, all=true)
	strname = string(name)
	(name in _symbols_excluded || startswith(strname, "_"))                                        && continue
	!(Base.isidentifier(name)  || (startswith(strname, "@") && Base.isidentifier(strname[2:end]))) && continue
	@eval export $name
end

end # module
