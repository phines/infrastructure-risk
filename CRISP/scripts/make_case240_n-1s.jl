using CSV;
include("..\\src\\ps_n-1_sec.jl")
include("..\\src\\CRISP_network.jl")

name = "data\\saved_ps\\case240_pserc_n-1";
case = "data\\case240_pserc\\"
ps = import_ps(case);
b1 = ps.branch[54,:];
append!(ps.branch,b1)
b2 = ps.branch[55,:];
append!(ps.branch,b2)
b3 = ps.branch[118,:];
append!(ps.branch,b3)
b4 = ps.branch[171,:];
append!(ps.branch,b4)
b5 = ps.branch[396:448,:];
append!(ps.branch,b5)
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
