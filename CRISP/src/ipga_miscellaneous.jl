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

module ipga_miscellaneous

using JuMP, LinearAlgebra, SparseArrays, DataFrames


@doc """`(optimal, feasible) = get_model_exit_status(M)`

Extract the post-optimized status of a JuMP Model `M`.

If a global/local optimum solution was found, then `optimal` is true.\n
If Model constraints were satisfied or nearly satisfied to the given tolerance, then `feasible` is true.""" ->
function get_model_exit_status(M::Model)
	termination_status_code = termination_status(M)
	result_status_code      = primal_status(M)

	optimal  = termination_status_code == MOI.OPTIMAL        || termination_status_code == MOI.LOCALLY_SOLVED
	feasible = result_status_code      == MOI.FEASIBLE_POINT || result_status_code      == MOI.NEARLY_FEASIBLE_POINT

	return optimal, feasible
end


@doc """`bit_vector = is_member(a, b)`

An implementation of MATLAB's `ismember` function!""" ->
function is_member(a, b)
	return .!isnothing.(indexin(a, b))
end


@doc """`i = find_in(a, b)`

An implementation of findin() from Julia versions <= 0.6.x.""" ->
function find_in(a, b)
	return findall(.!isnothing.(indexin(a, b)))
end


@doc """`(i, j) = index_in(a, b)`

A variant of the built-in function `indexin`, where `i` is the non-nothing version of the array returned by indexin,
and `j` is a bit-vector indicating true/false for each element in `a` present/absent in `b`.""" ->
function index_in(a, b)
	i = indexin(a, b)
	j = .!isnothing.(i)

	return i[j], j
end


@doc """`sparse_I = sp_eye(N, T)`

An implementation of speye() from Julia versions <= 0.6.x.

Input-\n
N: matrix size\n
T: the type of speye desired""" ->
function sp_eye(N::Int64; T::Type=Float64)
	return SparseMatrixCSC{T}(I,N,N)
end


@doc """`df_add_columns!(df, required_columns)`

Check the DataFrame `df` for the existence of `required_columns`, where the latter is a tuple of column names and data types, as illustrated below:

`(("id",    `Union{Missing,Int64}),
("kind",   `Union{Missing,String}),
("baseKV", `Union{Missing,Float64}),
("Vm",     `Union{Missing,Float64}),
("Va",     `Union{Missing,Float64}),
("Vsp",    `Union{Missing,Float64}))`

If a required column doesn't exist, then it is created.""" ->
function df_add_columns!(df::DataFrame, required_columns::Any)
	existing_columns = describe(df).variable

	N = nrow(df)

	for rc in required_columns
		if !in(Symbol(rc[1]), existing_columns)
			df[!,Symbol(rc[1])] = Array{rc[2]}(missing,N)
		elseif typeof(df[!,Symbol(rc[1])]) != rc[2] # then force typecast it
			df[!,Symbol(rc[1])] = Array{rc[2]}(df[!,Symbol(rc[1])])
		end
	end
end


@doc """`df_check_missing(df, column_names)`

Check the DataFrame `df` for missing data in the columns `column_names`. If they do exist, an error is thrown. """ ->
function df_check_missing(df::DataFrame, column_names::Array{String,1})
	n  = size(column_names,1)
	NA = falses(n)

	for i in 1:n
		if any(ismissing.(df[!,Symbol(column_names[i])]))
			NA[i] = true
			println(column_names[i])
		end
	end

	any(NA) && error("has missing data.")
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
