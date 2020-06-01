module TSSP

using SparseArrays
using LinearAlgebra

using JuMP

import SMPSReader
const SMPS = SMPSReader

"""

Reference formulation:
```
    min    c'x + q'y
    s.t.   A x       = b
           T x + W y = h
             x,    y ≥ 0
```

For the i-th scenario:
```
    min    c'x  + qi'y
    s.t.    A x        = b
           Ti x + Wi y = hi
              x,     y ≥ 0
```
where ``qi = q+ δq[i]``, ``hi = h + δh[i]``, ``Ti = T + ΔT[i]``, ``Wi = W + ΔW[i]``.
"""
struct TwoStageStochasticProgram
    m1::Int
    n1::Int
    m2::Int
    n2::Int

    # Template information
    A::SparseMatrixCSC{Float64,Int64}
    T::SparseMatrixCSC{Float64,Int64}
    W::SparseMatrixCSC{Float64,Int64}
    c::Vector{Float64}
    q::Vector{Float64}
    b::Vector{Float64}
    h::Vector{Float64}

    # Probabilistic data
    # This assumes a finite number of scenarios
    nscenarios::Int  # number of scenarios
    ΔTs::Vector{SparseMatrixCSC{Float64,Int64}}
    ΔWs::Vector{SparseMatrixCSC{Float64,Int64}}
    δqs::Vector{SparseVector{Float64}}
    δhs::Vector{SparseVector{Float64}}
    probs::Vector{Float64}
end

function all_realizations(X::SMPS.ScalarDiscrete)
    return [([(X.row_name, X.col_name, x)], p) for (x, p) in zip(X.support, X.p)]
end

function all_realizations(X::SMPS.BlockDiscrete)
    return zip(X.support, X.p)
end

"""
    TwoStageStochasticProgram(cdat, tdat, sdat)

Build a TwoStageStochasticProgram from smps data.
"""
function TwoStageStochasticProgram(cdat::SMPS.QPSReader.QPSData, tdat::SMPS.TimeData, sdat::SMPS.StocData)

    # Parition rows and columns into 1st and 2nd time periods
    j1 = cdat.varindices[tdat.cols[1]]  # Index of first 1st-period variable
    j2 = cdat.varindices[tdat.cols[2]]  # Index of first 2nd-period variable
    i1 = cdat.conindices[tdat.rows[1]]  # Index of first 1st-period constraint
    i2 = cdat.conindices[tdat.rows[2]]  # Index of first 2nd-period constraint
    @assert i1 == 1
    @assert j1 == 1
    @assert i2 <= cdat.ncon
    @assert j2 <= cdat.nvar

    con1 = i1:(i2 - 1)
    con2 = i2:cdat.ncon
    var1 = j1:(j2 - 1)
    var2 = j2:cdat.nvar
    m1, n1 = length(con1), length(var1)
    m2, n2 = length(con2), length(var2)

    # Convert to standard form
    # TODO: handle variable bounds as well
    nslack1 = 0
    nslack2 = 0
    srows = Int[]
    scols = Int[]
    svals = Float64[]
    rhs = zeros(m1 + m2)
    for (i, (l, u)) in enumerate(zip(cdat.lcon, cdat.ucon))
        first_period = i <= m1
        if l == -Inf && isfinite(u)
            # a'x ≤ u
            nslack1 += first_period
            nslack2 += !first_period
            push!(srows, i)
            push!(scols, nslack1 + nslack2)
            push!(svals, 1.0)
            rhs[i] = u
        elseif isfinite(l) && u == Inf
            # a'x ≥ l
            nslack1 += first_period
            nslack2 += !first_period
            push!(srows, i)
            push!(scols, nslack1 + nslack2)
            push!(svals, -1.0)
            rhs[i] = l
        elseif l == u
            # a'x = b
            rhs[i] = l
        else
            error("Unsupported bounds for row $i: [$l, $u]")
        end
    end

    M = sparse(cdat.arows, cdat.acols, cdat.avals, cdat.ncon, cdat.nvar)
    S = sparse(srows, scols, svals, m1+m2, nslack1+nslack2)

    # Build template data
    A = hcat(M[con1, var1], S[con1, 1:nslack1])
    T = hcat(M[con2, var1], spzeros(m2, nslack1))
    W = hcat(M[con2, var2], S[con2, (nslack1+1):end])
    c = zeros(n1 + nslack1); c[1:n1] .= cdat.c[1:n1]
    q = zeros(n2 + nslack2); q[1:n2] .= cdat.c[(n1+1):end]
    b = rhs[con1]
    h = rhs[con2]

    # Extract all scenarios
    R_indep = all_realizations.(sdat.indeps);
    R_blocks = all_realizations.(sdat.blocks);
    R_all = Base.Iterators.product(R_indep..., R_blocks...)
    nscenarios = length(R_all)

    # Construct perturbations
    ΔTs = SparseMatrixCSC{Float64,Int64}[]
    ΔWs = SparseMatrixCSC{Float64,Int64}[]
    δqs = SparseMatrixCSC{Float64,Int64}[]
    δhs = SparseMatrixCSC{Float64,Int64}[]
    probs = ones(nscenarios)
    for (k, r) in enumerate(R_all)
        pk = 1.0

        # We build the stochastic perturbation in COO format
        # sparse matrices/vectors are instantiated later
        trows = Int[]
        tcols = Int[]
        tvals = Float64[]

        wrows = Int[]
        wcols = Int[]
        wvals = Float64[]

        qind = Int[]
        qval = Float64[]

        hind = Int[]
        hval = Float64[]

        for (r_, p_) in r
            pk *= p_
            for (cname, vname, z) in r_
                if cdat.objname == cname
                    i = 0
                elseif haskey(cdat.conindices, cname)
                    i = cdat.conindices[cname]
                else
                    error("Unknown row $cname")
                end
    
                if cdat.rhsname == vname
                    j = 0
                elseif haskey(cdat.varindices, vname)
                    j = cdat.varindices[vname]
                else
                    error("Unknown variable $vname")
                end
    
                # Update correct coefficient
                if i == 0 && j > 0
                    # Objective
                    @assert j > n1
                    push!(qind, j-n1)
                    push!(qval, z - q[j - n1])
                elseif j == 0 && i > 0
                    # Right-hand side
                    @assert i > m1
                    push!(hind, i - m1)
                    push!(hval, z - h[i - m1])
                elseif i > 0 && j > 0
                    @assert i > m1
                    # Matrix coefficient
                    if j <= n1
                        push!(trows, i - m1)
                        push!(tcols, j)
                        push!(tvals, z - T[i - m1, j])
                    else
                        push!(wrows, i - m1)
                        push!(wcols, j - n1)
                        push!(wvals, z - W[i - m1, j - n1])
                    end
                else
                    error("Invalid coefficient pair ($i, $j)")
                end
            end
        end

        push!(ΔTs, sparse(trows, tcols, tvals, m2, n1+nslack1))
        push!(ΔWs, sparse(wrows, wcols, wvals, m2, n2+nslack2))
        push!(δqs, sparsevec(qind, qval, n2+nslack2))
        push!(δhs, sparsevec(hind, hval, m2))
        probs[k] = pk
    end

    # Done
    return TwoStageStochasticProgram(
        m1, n1 + nslack1, m2, n2 + nslack2,
        A, T, W, c, q, b, h,
        nscenarios, ΔTs, ΔWs, δqs, δhs, probs
        )
end

function deterministic_problem(tssp::TwoStageStochasticProgram, optimizer)
    model = JuMP.Model(optimizer)

    @variable(model, x[1:tssp.n1] >= 0)
    @variable(model, y[1:tssp.n2, 1:tssp.nscenarios] >= 0)
    @constraint(model, tssp.A * x .== tssp.b)
    q = zeros(tssp.n2, tssp.nscenarios)  # This will be the objective vector for y

    # Second period constraints
    for (k, (ΔT, ΔW, δq, δh, p)) in enumerate(zip(tssp.ΔTs, tssp.ΔWs, tssp.δqs, tssp.δhs, tssp.probs))
        @constraint(model, (tssp.T + ΔT) * x + (tssp.W + ΔW) * y[:, k] .== (tssp.h + δh))
        # Compute correct objective
        q[:, k] .= (p .* (tssp.q .+ δq))
    end

    # Set the objective
    @objective(model, Min, dot(tssp.c, x) + dot(q[:], y[:]))

    return model
end

end # module
