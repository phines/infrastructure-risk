using CSV;
include("..\\src\\ps_n-1_s_gen.jl")
include("..\\src\\CRISP_network_gen.jl")

name = "data\\saved_ps\\case73_n-1";
case = "data\\case73"
name1 = "data\\saved_ps\\case73\\"
if isdir(name1)
else
    mkdir(name1)
end
ps = import_ps(case);
ps.gen = ps.gen[1:end-1,:]
export_ps(ps,name1)
ps.storage = DataFrame(bus = [313.0], E = [75.0], Ps = [0.0], Emax = [150.0], Emin = [0.0], Psmax = [50.0], Psmin = [-50.0], status = [1]);
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
