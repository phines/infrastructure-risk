using SparseArrays
using LinearAlgebra

A = sprand(10,10,.9)
b_dense = Array(sprand(10,1,.8)) # works
b_sparse = sprand(10,1,.8) # doesn't work
x = A\b_dense # works
x = A\b_sparse # error
lu!(A) # error

