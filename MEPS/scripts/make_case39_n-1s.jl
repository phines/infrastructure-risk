using CSV;
include("..\\src\\ps_n-1_sec.jl")
include("..\\src\\CRISP_network.jl")

name = "data\\saved_ps\\case39_n-1";
case = "data\\case39\\"
ps = import_ps(case);
b1 = ps.branch[5,:];
append!(ps.branch,b1)
b2 = ps.branch[20,:];
append!(ps.branch,b2)
b3 = ps.branch[27,:];
append!(ps.branch,b3)
b4 = ps.branch[33,:];
append!(ps.branch,b4)
b5 = ps.branch[34,:];
append!(ps.branch,b5)
b6 = ps.branch[37,:];
append!(ps.branch,b6)
b7 = ps.branch[39,:];
append!(ps.branch,b7)
b8 = ps.branch[41,:];
append!(ps.branch,b8)
b9 = ps.branch[46,:];
append!(ps.branch,b9)
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
