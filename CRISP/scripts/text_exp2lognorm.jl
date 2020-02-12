using Random
using DataFrames
using CSV
rng = MersenneTwister(10)
include("../src/CRISP_initiate.jl")
N = 100000
mu = 3.66
sigma = 2.43
P1 = 0.1:0.1:0.5
p2 = 10
factor = 1.1
rest_times1 = RecoveryTimes(mu,sigma,N)
for p1 in P1
    rest_times2 = exp2lognorm(N,p1,p2,factor)
    results = DataFrame(log_normal_dist = rest_times1, weibull_gone_bad = rest_times2)
    CSV.write("results\\test_exponential2lognormal_p1_$p1.csv",results)
end
p1 = 0.3
P2 = 15:5:45
for p2 in P2
    rest_times2 = exp2lognorm(N,p1,p2,factor)
    results = DataFrame(log_normal_dist = rest_times1, weibull_gone_bad = rest_times2)
    CSV.write("results\\test_exponential2lognormal_p2_$p2.csv",results)
end
