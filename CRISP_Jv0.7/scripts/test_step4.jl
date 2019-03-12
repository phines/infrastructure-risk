#using CRISP_LSOPF
include("..\\src\\CRISP_RT.jl")
Restore = include("test_step1+step2+step3.jl") # data frame [times, load shed in cost per hour]
ResilienceTri = crisp_res(Restore);
