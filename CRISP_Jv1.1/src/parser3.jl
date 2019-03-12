using DataFrames
using CSV

function mp2ps(case, ybus_includeZshunts = false)
    (mpBusData, numBusRows, nunBusCols) = readBusData(case);
    (mpBranchData, numBranchRows, numBranchCols)= readBranchData(case);
    (mpGenData, numGenRows, numGenCols) = readGenData(case)
    (mpGenCostData, numGenCostRows, numGenCostCols) = readGenCostData(case);
    mpBaseMVA = readBaseMVA(case);

    psGenCostData = makeGenCostDF();
    shuntdf = makeShuntDF();

    psBusData = mpBus2psBus(mpBusData);
    psBranchData = mpBranch2psBranch(mpBranchData);
    psGenData = mpGen2psGen(mpGenData);
    psShuntData = pushToDataFrame(makeShuntDF(), getShuntData(mpBusData));

    psCase = PSCase(mpBaseMVA, psBusData, psBranchData, psGenData, psGenCostData, psShuntData);

    return psCase
end

function findBusData(filename)
    busDataStarts = 0;
    busDataEnds = 0;
    inBusData = false;
    open(filename) do file
        for i in enumerate(eachline(file))
            if (occursin(r".*\.bus\s*=\s*\[", i[2]))
                busDataStarts = i[1];
                inBusData = true;
            end
            if (occursin(r".*\](;)?", i[2]) & inBusData)
                busDataEnds = i[1];
                inBusData = false;
                break;
            end
        end
    end
    return busDataStarts, busDataEnds;
end

function findBranchData(filename)
    branchDataStarts = 0;
    branchDataEnds = 0;
    inBranchData = false;
    open(filename) do file
        for i in enumerate(eachline(file))
            if (occursin(r".*\.branch\s*=\s*\[", i[2]))
                branchDataStarts = i[1];
                inBranchData = true;
            end
            if (occursin(r".*\](;)?", i[2]) & inBranchData)
                branchDataEnds = i[1];
                inBranchData = false;
                break;
            end
        end
    end
    return branchDataStarts, branchDataEnds;
end

function findGenCostData(filename)
    genCostDataStarts = 0;
    genCostDataEnds = 0;
    inGenCostData = false;
    open(filename) do file
        for i in enumerate(eachline(file))
            if (occursin(r".*\.gencost\s*=\s*\[", i[2]))
                genCostDataStarts = i[1];
                inGenCostData = true;
            end
            if (occursin(r".*\](;)?", i[2]) & inGenCostData)
                genCostDataEnds = i[1];
                inGenCostData = false;
                break;
            end
        end
    end
    return genCostDataStarts, genCostDataEnds;
end

function findGenData(filename)
    genDataStarts = 0;
    genDataEnds = 0;
    inGenData = false;
    open(filename) do file
        for i in enumerate(eachline(file))
            if (occursin(r".*\.gen\s*=\s*\[", i[2]))
                genDataStarts = i[1];
                inGenData = true;
            end
            if (occursin(r".*\](;)?", i[2]) & inGenData)
                genDataEnds = i[1];
                inGenData = false;
                break;
            end
        end
    end
    return genDataStarts, genDataEnds;
end

function findBaseMVA(filename)
    baseMVAisAt = 0;
    open(filename) do file
        for i in enumerate(eachline(file))
            if (occursin(r".*\.baseMVA\s*=\s*", i[2]))
                baseMVAisAt = i[1];
                break;
            end
        end
    end
    return baseMVAisAt;
end

function findAreaNames(filename)
    #TODO
end

function findBusNames(filename)
    #TODO
end



function readBusData(filename)
    # fails if a line has only a comment.
    busData = [];
    (busDataStart, busDataEnd) = findBusData(filename);
    numRows = busDataEnd - busDataStart - 1;
    numCols = 13;
    open(filename) do file
        for i in enumerate(eachline(file))
            if (i[1]> busDataStart && i[1] < busDataEnd)
                line = i[2];
                if (occursin(r"%.*", line))
                    line = replace(line, r"%.*", " ");
                end
                row = map(x -> parse(ComplexF64,x), collect(m.match for m in eachmatch(r"((?:[-+]?[0-9]*(e-|\.)?[0-9]+)+)", line)));
                for data in row
                    push!(busData, data);
                end
            end
        end
    end
    busData = transpose(reshape(busData, (numCols, numRows)));
    return busData, numRows, numCols;
end

function readBranchData(filename)
    # fails if a line has only a comment.
    branchData = [];
    (branchDataStart, branchDataEnd) = findBranchData(filename);
    numRows = branchDataEnd - branchDataStart - 1;
    numCols = 13;
    open(filename) do file
        for i in enumerate(eachline(file))
            if (i[1]> branchDataStart && i[1] < branchDataEnd)
                line = i[2];
                if (occursin(r"%.*", line))
                    line = replace(line, r"%.*", " ");
                end
                row = map(parse, collect(m.match for m in eachmatch(r"((?:[-+]?[0-9]*(e-|\.)?[0-9]+)+)", line)));
                for data in row
                    push!(branchData, data);
                end
            end
        end
    end
    branchData = transpose(reshape(branchData, (numCols, numRows)));
    return branchData, numRows, numCols;
end

function readGenData(filename)
      # fails if a line has only a comment.
    genData = [];
    (genDataStart, genDataEnd) = findGenData(filename);
    numRows = genDataEnd - genDataStart - 1;
    numCols = 21;
    open(filename) do file
        for i in enumerate(eachline(file))
            if (i[1]> genDataStart && i[1] < genDataEnd)
                line = i[2];
                if (occursin(r"%.*", line))
                    line = replace(line, r"%.*", " ");
                end
                row = map(parse, collect(m.match for m in eachmatch(r"((?:[-+]?[0-9]*(e-|\.)?[0-9]+)+)", line)));
                for data in row
                    push!(genData, data);
                end
            end
        end
    end
    genData = transpose(reshape(genData, (numCols, numRows)));
    return genData, numRows, numCols;
end

function readGenCostData(filename)
      # fails if a line has only a comment.
    genCostData = [];
    (genCostDataStart, genCostDataEnd) = findGenCostData(filename);
    numRows = genCostDataEnd - genCostDataStart - 1;
    numCols = 0;
    open(filename) do file
        for i in enumerate(eachline(file))
            if (i[1]> genCostDataStart && i[1] < genCostDataEnd)
                line = i[2];
                if (occursin(r"%.*", line))
                    line = replace(line, r"%.*", " ");
                end
                row = map(parse, collect(m.match for m in eachmatch(r"((?:[-+]?[0-9]*(e-|\.)?[0-9]+)+)", line)));
                for data in row
                    push!(genCostData, data);
                end
            end
        end
    end
    numCols = convert(Int64,length(genCostData)/numRows);
    genCostData = transpose(reshape(genCostData, (numCols, numRows)));
    return genCostData, numRows, numCols;
end

function readBaseMVA(filename)
    baseMVA = 1;
    baseMVAisAt = findBaseMVA(filename)
    open(filename) do file
        for i in enumerate(eachline(file))
            if (i[1] == baseMVAisAt)
                line = i[2];
                baseMVAString = collect(m.match for m in eachmatch(r".*\.baseMVA\s*=\s*((?:[-+]?[0-9]*(e-|\.)?[0-9]+)+);", line));
                baseMVA = parse(baseMVAString.captures[1]);
                break;
            end
        end
    end
    return baseMVA
end

function getShuntData(busData, ybus_includeZshunts = false)
    #Find non zero Pd + j*Qd buses from busData;
        #return values and indices; assign to S and ixS respectively
    #Find non zero Gs and Ps buses from busData;
        #return values and indices; assign to Y and ixY respectively
    # ixY + ixS = numOfRows for Shunt DataType
    # numOfCols = num of fields in Shunt DF
    # intiliaze shuntData with zeros(numOfRows, numOfCols)
    # assign values to fields
    # return shunt data

    ##############################
    # PS Constants
    fixed = 1;
    switched = 0;
    ISO = 4;
    ##############################

    Pd = busData[:, 3];
    Qd = busData[:, 4];
    Gs = busData[:, 5];
    Bs = busData[:, 6];

    ixS = findall((Pd + im*Qd) .!= 0);
    ixY = findall((Gs + im*Bs) .!= 0);

    S = Pd[ixS] + im*Qd[ixS];
    Y = Gs[ixY] + im*Bs[ixY];

    ns = length(ixS);
    ny = length(ixY);
    nsy = ns + ny;
    numRows = nsy;
    numCols = 13;

    shuntData = zeros(numRows, numCols);

    shuntData[:, 1] .= busData[[ixS;ixY], 1]; # Bus ID
    shuntData[:, 2] .= fixed; # Shunt Type: TODO
    shuntData[findall(busData[[ixS;ixY], 2] .!= 4), 3] .= 1; # Bus Type != 4 == Isolated Bus
    shuntData[(1:ns), 4] .= real(S);
    shuntData[(ns+1 : nsy) , 4] .= real(Y);
    shuntData[(1:ns), 5] .= imag(S);

    if (ybus_includeZshunts)
        shuntData[(ns+1 : nsy), 5] = imag(Y);
    else
        shuntData[(ns+1 : nsy), 5] = -imag(Y);
    end

    shuntData[(1:ns), 6] .= 1;
    shuntData[(ns+1 : nsy), 7] .= 1;

    return shuntData;
end


function mpBus2psBus(mpBusData)
    psBusData = makeBusDF();
    if (size(psBusData) != size(mpBusData))
        numOfCols = size(psBusData, 2) - size(mpBusData, 2);
        mpBusData = addColumns!(mpBusData, numOfCols)
    end
    psBusData = pushToDataFrame(psBusData, mpBusData);
    return psBusData;
end

function mpBranch2psBranch(mpBranchData)
    psBranchData = makeBranchDF();
    if (size(psBranchData) != size(mpBranchData))
        numOfCols = size(psBranchData, 2) - size(mpBranchData, 2);
        mpBranchData = addColumns!(mpBranchData, numOfCols)
    end
    psBranchData = pushToDataFrame(psBranchData, mpBranchData);
    return psBranchData;
end

function mpGen2psGen(mpGenData)
    psGenData = makeGenDF();
    if (size(psGenData) != size(mpGenData))
        numOfCols = size(psGenData, 2) - size(mpGenData, 2);
        mpGenData = addColumns!(mpGenData, numOfCols)
    end
    psGenData = pushToDataFrame(psGenData, mpGenData);
    return psGenData;
end

function mpGenCost2psGenCost(mpGenCostData)
    #TODO
end

mutable struct PSCase
    baseMVA::Int64
    bus::DataFrame
    branch::DataFrame
    gen::DataFrame
    gencost::DataFrame
    shunt::DataFrame
end

function makeBusDF()
    busdf = DataFrame(id = Int[],
            bus_type = Int[],
            Pd = Float64[],
            Qd = Float64[],
            Gs = Float64[],
            Bs = Float64[],
            area = Float64[],
            Vm = Float64[],
            Va = Float64[],
            baseKV = Float64[],
            zone = Float64[],
            Vmax = Float64[],
            Vmin = Float64[],
            LAM_P = Float64[],
            LAM_Q = Float64[],
            MU_VMAX = Float64[],
            MU_VMIN = 1Float64[],
            locX = Float64[],
            locY = Float64[],)
    return busdf
end

function makeBranchDF()
    branchdf = DataFrame(f = Int[],
            t = Int[],
            R = Float64[],
            X = Float64[],
            B = Float64[],
            rateA = Float64[],
            rateB = Float64[],
            rateC = Float64[],
            tap = Float64[],
            shift = Float64[],
            status = Float64[],
            angMin = Float64[],
            angMax = Float64[],
            Pf = Float64[],
            Qf = Float64[],
            Pt = Float64[],
            Qt = Float64[],
            MU_SF = Float64[],
            MU_ST = Float64[],
            MU_ANGMIN = Float64[],
            MU_ANGMAX = Float64[],
            # ---these columns (22-27) are exclusive to our ps---
            ImF = Float64[],
            ImT = Float64[],
            switchable = Float64[],
            failProb = Float64[],
            branch_type = Float64[],
            id = Float64[])
    return branchdf
end

function makeGenDF()
    gendf = DataFrame(bus = Int[], #GEN_BUS
            Pg = Float64[],
            Qg = Float64[],
            Qmax = Float64[],
            Qmin = Float64[],
            Vsp = Float64[], #Vg
            mBase = Float64[],
            status = Float64[],
            Pmax = Float64[],
            Pmin = Float64[],
            Pc1 = Float64[],
            Pc2 = Float64[],
            Qc1min = Float64[],
            Qc1max = Float64[],
            Qc2min = Float64[],
            Qc2max = Float64[],
            rampAGC = Float64[], #RAMP_AGC
            ramp10 = Float64[], #RAMP_10
            ramp30 = Float64[], #RAMP_30
            rampQ = Float64[], #RAMP_Q
            npf = Float64[], #APF
            MU_PMAX = Float64[],
            MU_PMIN = Float64[],
            MU_QMAX = Float64[],
            MU_QMIN = Float64[],
            #---these columns (26-27) are exclusive to our ps---#
            bus_type = Float64[], #bus-type (i.e., PV/PQ/PF/REF)
            cost = Float64[] #marginal-cost ($/MW) of power
            )
end

function makeShuntDF()
    shuntdf = DataFrame(
        bus = Int[], #BUS_I of the parent bus
        shuntType = Int[], #fixed/switched?
        status = Int[], #ON/OFF?
        #fixed shunts: e.g.,  loads and capacitors --
        P = Float64[], #P-demand (MW)   at 1 p.u. V
        Q = Float64[], #Q-demand (MVAr) at 1 p.u. V
        frac_S = Float64[], #fraction (of P+jQ) that's constant-power
        frac_Y = Float64[], #                          constant-impedance
        frac_E = Float64[], #                          an expontential-load (P*Vm^gamma)
        gamma = Float64[], #gamma -- see above
        #switched shunts: e.g., synchronous condensers and static VAr compensators (SVCs) --
        Qmin = Float64[], #device Q-minimum
        Qmax = Float64[], #       Q-maximum
        Qinc = Float64[], #       Q-increment
        Qcon = Float64[])#       Q-consumption
end

function makeGenCostDF()
    gencostdf = DataFrame(
        model = Int[],
        startup_cost = Float64[],
        shutdown_cost = Float64[],
        num_model_params = Float64[]
    )
end

function addColumns!(mat, numOfCols)
    numOfRows = size(mat, 1);
    addCols = zeros(numOfRows, numOfCols);
    bigMat = hcat(mat, addCols)
    return bigMat;
end

function pushToDataFrame(df, arr)
    numRows = size(arr, 1);
    for i in (1:numRows)
        df = push!(df, arr[i, :])
    end
    return df;
end
