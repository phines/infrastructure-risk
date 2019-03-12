using JuMP
using Gurobi
#using Ipopt

include("../src/parser.jl")

# load the case data
ps = mp2ps("C:\\Users\\mkellygo\\Documents\\Github\\infrastructure-risk\\CRISP\\data\\case6ww.m") #("../data/case6ww.m")
#ps = mp2ps("../data/case6ww.m")
ps.branch[2,:status]=0;
ps.branch[5,:status]=0;
function crisp_soc_lsopf(ps)
    ### collect the data that we will need ###
    # bus data
    n = size(ps.bus,1) # the number of buses
    bi = sparse(ps.bus[:id],fill(1,n),collect(1:n)) # helps us to find things
    vup = ps.bus[:,:Vmax];#upper voltage bounds
    vlo = ps.bus[:,:Vmin];#lower voltage bounds
    # load data
    nd = size(ps.shunt,1)
    D = bi[ps.shunt[:bus]]
    Pd = ps.shunt[:P] ./ ps.baseMVA .* ps.shunt[:status]
    Qd = ps.shunt[:Q] ./ ps.baseMVA .* ps.shunt[:status]
    # gen data
    ng = size(ps.gen,1)
    G = bi[ps.gen[:bus]]
    Pg = ps.gen[:Pg] ./ ps.baseMVA .* ps.gen[:status]
    Pmax = ps.gen[:Pmax] ./ ps.baseMVA .* ps.gen[:status]
    Pmin = ps.gen[:Pmin] ./ ps.baseMVA .* ps.gen[:status]
    Qmax = ps.gen[:Qmax] ./ ps.baseMVA .* ps.gen[:status]
    Qmin = ps.gen[:Qmin] ./ ps.baseMVA .* ps.gen[:status]
    # branch data
    brst = (ps.branch[:status].==1)
    F = bi[ps.branch[brst,:f]]
    T = bi[ps.branch[brst,:t]]
    nl = length(T);
    flow0 = ps.branch[brst,:Pf]./ps.baseMVA
    flow_max = ps.branch[brst,:rateA]./ps.baseMVA # this could also be rateB
    Xinv = (1 ./ ps.branch[brst,:X])
    Zinv = (1 ./ (ps.branch[brst,:R] + im.*ps.branch[brst,:X]))
    B = sparse(F,T,-Xinv,n,n) +
        sparse(T,F,-Xinv,n,n) +
        sparse(T,T,+Xinv,n,n) +
        sparse(F,F,+Xinv,n,n)
    bc = +ps.branch[brst,:B];
    # Admittance matrix
    Yij = sparse(F,T,-Zinv, n,n) +
          sparse(T,F,-Zinv, n,n) +
          sparse(T,T,+Zinv,n,n) +
          sparse(F,F,+Zinv,n,n)
    bij = real(-Zinv);
    gij = imag(-Zinv);
    Flows = sparse(F,T,1, n,n) +
        sparse(T,F,1, n,n);
    Ys = ps.bus[:,:Gs] + im.*ps.bus[:,:Bs]# line charging admittance matrix.
    ### Build the optimization model ###
    m = Model(solver = GurobiSolver())
    #m = Model(solver = IpoptSolver())
    # variables
    @variable(m,Pg[1:ng])
    @variable(m,Qg[1:ng])
    @variable(m,dTheta[1:n])
    @variable(m,Pij[1:nl])
    @variable(m,Qij[1:nl])
    @variable(m,Pji[1:nl])
    @variable(m,Qji[1:nl])
    @variable(m,W_ii[1:n])
    @variable(m,Wr_ij[1:nl])
    @variable(m,Wi_ij[1:nl])
    @variable(m,Ws_ii[1:n])
    @variable(m, z_v[1:n])
    @variable(m, z_g[1:ng])
    @variable(m, z_d[1:nd])
    @variable(m, z_s[1:n])
    # variable bounds
    @constraint(m, Pg .<= Pmax)
    @constraint(m, Pg .>= 0)
    @constraint(m, Qg .<= Qmax)
    @constraint(m, Qg .>= Qmin)
    @constraint(m, 0 .<= z_v)
    @constraint(m, 0 .<= z_g)
    @constraint(m, 0 .<= z_d)
    @constraint(m, 0 .<= z_s)
    @constraint(m, 1 .>= z_v)
    @constraint(m, 1 .>= z_g)
    @constraint(m, 1 .>= z_d)
    @constraint(m, 1 .>= z_s)
#    @constraint(m,dTheta[1] == 0)
    @constraint(m, constr[i=1:n], W_ii[i] <= z_v[i]*vup[i]^2)
    @constraint(m, constr[i=1:n], W_ii[i] >= z_v[i]*vlo[i]^2)
    #@constraint(m, (Wr_ij.^2+Wi_ij.^2) .<= W_ii[F].*W_ii[T])
    # McCormick envelope representing Ws_ii = z_s.*Wii
    #@constraint(m, constr[i=1:n], Ws_ii[i] >= z_s[i]*vlo[i]^2)
    #@constraint(m, constr[i=1:n], Ws_ii[i] >= W_ii[i] + z_s[i]*vup[i]^2 - vup[i]^2)
    #@constraint(m, constr[i=1:n], Ws_ii[i] <= z_s[i]*vup[i]^2)
    #@constraint(m, constr[i=1:n], Ws_ii[i] <= W_ii[i] + z_s[i]*vlo[i]^2 - vlo[i]^2)
    # objective: serve as much load as possible and keep as much of the network connected as possible
    C = 100;
    @objective(m,Min,-sum(z_v)-sum(z_g)-sum(z_s)-C*sum(z_d.*Pd))
    # mapping matrix to map loads/gens to buses
    M_D = sparse(D,1:nd,1.0,n,nd)
    M_G = sparse(G,1:ng,1.0,n,ng)
    # Power balance equality constraint
    Ysr = real(Ys);
    Ysi = imag(Ys);
    @constraint(m, constr[i=1:n], (M_G*Pg)[i] - (M_D*z_d)[i].*(M_D*Pd)[i]
                     - Ysr[i].*W_ii[i] ==  sum(Pij[F.==i]) - sum(Pji[T.==i]))
    @constraint(m, constr[i=1:n], (M_G*Qg)[i] - (M_D*z_d)[i].*(M_D*Qd)[i]
                    - Ysi[i].*W_ii[i] ==  sum(Qij[F.==i]) - sum(Qji[T.==i])) #
    # Power flow constraints
    @constraint(m, Pij.^2+Qij.^2 .<= flow_max.^2)
    @constraint(m, Pji.^2+Qji.^2 .<= flow_max.^2)
    @constraint(m, Pij .== (bij).*W_ii[F] - (bij.*Wr_ij+gij.*Wi_ij))
    @constraint(m, Pji .== (bij).*W_ii[T] + (bij.*Wr_ij-gij.*Wi_ij))
    @constraint(m, Qij .== (-gij - bc/2).*W_ii[F] + (gij.*Wr_ij-bij.*Wi_ij))
    @constraint(m, Qji .== (-gij - bc/2).*W_ii[T] + (gij.*Wr_ij+bij.*Wi_ij))
    @constraint(m,tan(-pi/2).*Wr_ij .<= Wi_ij) 
    @constraint(m,tan(pi/2).*Wr_ij .>= Wi_ij)
    ### solve the model ###
    solve(m)
    # collect/return the outputs
    dPd_star = getvalue(z_d).*Pd.*ps.baseMVA
    dPg_star = getvalue(z_g).*Pg.*ps.baseMVA
    return (dPd_star, dPg_star)
end
