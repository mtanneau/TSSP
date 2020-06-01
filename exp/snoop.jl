include("colgen.jl")

const TSSP_DIR = joinpath(@__DIR__, "../data/tssp")

tssp = read_tssp(
    joinpath(TSSP_DIR, "4node.cor"),
    joinpath(TSSP_DIR, "4node.tim"),
    joinpath(TSSP_DIR, "4node.sto.128"),
)

run_colgen(tssp, SAVE_RMP=false)