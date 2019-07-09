using CSV;
include("..\\src\\ps_n-1_s_gen.jl")
include("..\\src\\CRISP_network_gen.jl")

name = "data\\saved_ps\\case73_n-1_noPV_noWind_noStorage";
case = "data\\case73\\"
ps = import_ps(case);
ps.gen = ps.gen[1:96,:]
b1 = ps.branch[52,:];
append!(ps.branch,b1)
b2 = ps.branch[90,:];
append!(ps.branch,b2)
secure_ps!(ps);
# save ps structure
if isdir(name)
else
    mkdir(name)
end
export_ps(ps,name)
#save case info
line1 = "$case made n-1 secure";
write("$name\\case_info.txt",line1);
