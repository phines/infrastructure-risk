# Step 4 of Crisp model. Measure of resilience triangle.
using DataFrames
using LinearAlgebra

function crisp_res(Restore)
  times = Restore[1] = times # time at which restoration occurs
  cost_load_shed = Restore[2]; # cost per hour
  T = times[2:end]-times[1:end-1]
  res = sum(T.*cost_load_shed[1:end-1])
  return res
end
