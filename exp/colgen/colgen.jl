using LinearAlgebra
using Printf
using SparseArrays
using Random

import Linda
import TSSP

import Gurobi
const GRBENV = Gurobi.Env()

import SMPSReader
const SMPS = SMPSReader

using JuMP

const RMPDIR = joinpath(@__DIR__, "../../data/rmp")

"""
    LPOracle

Oracle for LP sub-problems.
"""
mutable struct LPOracle <: Linda.Oracle
    index::Int
    
    m2::Int
    n2::Int
    
    sp::JuMP.Model
    h::Vector{Float64}
    T::SparseMatrixCSC{Float64,Int64}
    ΔT::SparseMatrixCSC{Float64,Int64}
    
    # Shadow prices
    farkas::Bool
    π::Vector{Float64}
    σ::Float64
    
    function LPOracle(idx::Int, m2::Int, n2::Int, T, ΔT, W, ΔW, q, h, p::Float64, optimizer)
        sp = JuMP.direct_model(MOI.instantiate(optimizer))
        
        @variable(sp, θ[1:m2])
        @constraint(sp, W'θ .+ ΔW'θ .<= (p .* q))
        
        @objective(sp, Min, h'θ)
        
        return new(idx, m2, n2, sp, h, T, ΔT, false, zeros(size(T, 2)), 0.0)
    end
end

import Linda:update!, optimize!, get_columns, get_dual_bound

function Linda.update!(oracle::LPOracle, farkas::Bool, π::Vector{Float64}, σ::Float64)

    # Update shadow prices
    oracle.farkas = farkas
    oracle.π .= π
    oracle.σ  = σ
    
    # Update objective
    θ = oracle.sp[:θ]
    @objective(oracle.sp, Min, dot( (!farkas .* oracle.h) - oracle.T * π - oracle.ΔT * π, θ))
        
    return nothing
end

function Linda.optimize!(oracle::LPOracle)
    JuMP.optimize!(oracle.sp)
    return nothing
end

function Linda.get_columns(oracle::LPOracle)
    
    # TODO: check solution status
    θ = oracle.sp[:θ]
    st = JuMP.termination_status(oracle.sp)
    
    st == JuMP.MOI.OPTIMAL || st == JuMP.MOI.DUAL_INFEASIBLE || error("Sub-problem exited with status $st")
    JuMP.has_values(oracle.sp) || error("SP solve with status $st but no values available")
    
    is_vertex = (st == JuMP.MOI.OPTIMAL )
    θ_ = value.(θ)
    
    # Compute real cost
    cost = dot(oracle.h, θ_)
    
    # Get reduced cost
    rc = is_vertex ? (JuMP.objective_value(oracle.sp) - oracle.σ) : -Inf
    
    u = sparsevec(oracle.T' * θ_ + oracle.ΔT' * θ_)
    
    
    # Build column
    col = Linda.Column(
        cost,
        oracle.index,
        u.nzind, u.nzval,
        is_vertex
    )
    
    return [(col, rc)]
end

function get_dual_bound(oracle::LPOracle)
    if oracle.farkas
        return -Inf
    else
        st = JuMP.termination_status(oracle.sp)
        if st == JuMP.MOI.OPTIMAL
            return JuMP.objective_value(oracle.sp)
        elseif st == JuMP.MOI.DUAL_INFEASIBLE
            return JuMP.objective_value(oracle.sp)
        else
            error("Cannot query dual bound: SP exited with status $st")
        end
    end
end


function read_tssp(fcore::String, ftime::String, fstoc::String)
    cdat = SMPS.read_core_file(fcore)
    tdat = SMPS.read_time_file(ftime)
    sdat = SMPS.read_stoch_file(fstoc)

    return TSSP.TwoStageStochasticProgram(cdat, tdat, sdat)
end

"""
    generate_oracles(tssp)

Generate one LPOracle for each scenario.
"""
generate_oracles(tssp) = [
    LPOracle(
        k, tssp.m2, tssp.n2,
        tssp.T, tssp.ΔTs[k],
        tssp.W, tssp.ΔWs[k],
        tssp.q + tssp.δqs[k],
        - (tssp.h + tssp.δhs[k]),
        tssp.probs[k],
        () -> Gurobi.Optimizer(GRBENV)
    )
    for k in 1:tssp.nscenarios
]

"""
    generate_master(tssp)

Generate the initial master problem.


"""
function generate_master(tssp, oracles; M=1e4)
    rmp = Gurobi.Optimizer(GRBENV, OutputFlag=0, Threads=1)
    jrmp = JuMP.direct_model(rmp)
    
    # Instantiate RMP
    @variable(jrmp, ηp[1:tssp.m1] >= 0)
    @variable(jrmp, ηm[1:tssp.m1] >= 0)
    @variable(jrmp,  ξ[1:tssp.n1] >= 0)
    @variable(jrmp,  ζ[1:tssp.n1] >= 0)
    
    # Linking constraints
    @constraint(jrmp, clink, (tssp.A' * ηp) .- (tssp.A' * ηm) .+ ξ .- ζ .== tssp.c)
    
    # Convexity constraints
    @constraint(jrmp, ccvx, 0 .== ones(tssp.nscenarios))
    
    # Objective
    @objective(jrmp, Min, -dot(tssp.b, ηp) + dot(tssp.b, ηm) + M*sum(ζ))
    
    mp = Linda.Master(
        rmp, getproperty.(ccvx, :index), getproperty.(clink, :index), tssp.c,
        vcat(getproperty.(ηp, :index), getproperty.(ηm, :index), getproperty.(ξ, :index), getproperty.(ζ, :index)),
        MOI.VariableIndex[]
    )

    # Add initial columns
    π0 = zeros(tssp.n1)
    σ0 = 0.0
    cols = Vector{Linda.Column}(undef, length(oracles))

    for (k, oracle) in enumerate(oracles)
        JuMP.set_optimizer_attribute(oracle.sp, "OutputFlag", 0)
        JuMP.set_optimizer_attribute(oracle.sp, "Threads", 1)
        
        Linda.update!(oracle, false, π0, σ0)
        # We minimize 0 to ensure we get a vertex
        @objective(oracle.sp, Min, 0)
        Linda.optimize!(oracle)
        col, rc = Linda.get_columns(oracle)[1]
        
        cols[k] = col
    end

    for col in cols
        Linda.add_column!(mp, col)
    end
    
    return mp
end

function run_colgen(tssp; fname="RMP", SAVE_RMP::Bool=false)
    env = Linda.LindaEnv()
    env.verbose.val = true

    # Generate at most |S| / 10 columns at each outer iteration
    env.num_columns_max.val = div(tssp.nscenarios, 5)

    # Generate oracles
    oracles = generate_oracles(tssp)
    mp = generate_master(tssp, oracles)

    cglog = Dict()

    # Callback function to save the current RMP
    # every 10 iterations
    function cg_callback()
        if SAVE_RMP
        niter::Int = cglog[:n_cg_iter]
            if niter > 0 && niter %10 == 0
                frmp = joinpath(RMPDIR, "$(fname)_$(tssp.nscenarios)_$(niter).mps.bz2")
                @info "Saving RMP at $frmp"
                # TODO: file name
                Gurobi.write_model(mp.rmp.inner, frmp)
            end
        end
        return nothing
    end

    # Solve the problem by column generation
    Random.seed!(42)
    @time Linda.solve_colgen!(env, mp, oracles, cg_log=cglog, callback=cg_callback)

    # Save the last RMP
    if SAVE_RMP
        frmp = joinpath(RMPDIR, "$(fname)_$(tssp.nscenarios)_$(cglog[:n_cg_iter]).mps.bz2")
        @info "Saving last RMP at $frmp"
        Gurobi.write_model(mp.rmp.inner, frmp)
    end

    # Done
    @info "Colgen stats" mp.primal_lp_bound mp.dual_bound
    return nothing
end

function main(args)
    tssp = read_tssp(args[1], args[2], args[3])

    fname = length(args) >= 4 ? args[4] : "RMP"
    run_colgen(tssp, fname=fname, SAVE_RMP=true) 

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end