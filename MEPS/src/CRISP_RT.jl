# Step 4 of Crisp model. Measure of resilience triangle.
using DataFrames
using LinearAlgebra

function crisp_res(Restore)
  times = Restore[1]; # time at which restoration occurs
  cost_load_shed = Restore[2]; # cost per hour
  T = times[3:end]-times[2:end-1];
  res = sum(T.*cost_load_shed[2:end-1])/60; #converting to hours from minutes.
  return res
end
