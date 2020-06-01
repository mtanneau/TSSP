# TSSP

Data and scripts for generating a collection of master problems from Two-Stage Stochastic Programming instances.

## Installation

1. Clone/download this repository
1. Download and install Gurobi
1. Install packages
    ```bash
    julia --project -e 'using Pkg; Pkg.instantiate()'
    ```

## Instances

We consider the TSSP instances with at least 1000 scenarios, which are displayed below.

| `.cor` file | `.tim` file | `.sto` file |
|:------------|:------------|:------------|
`env.cor`, `env.cor.diss` | `env.tim` | `env.sto.N`
`4node.cor`, `4node.cor.base` | `4node.tim` | `4node.sto.N`
`phone.cor` | `phone.tim` | `phone.sto`
`stormG2.cor` | `stormG2.tim` | `stormG2_1000.sto`

We remove `assets` and `phone` instances because they are solved in a single iteration.

Each instance is solved by column-generation, with Gurobi (default parameters) as the MP solver.
Sub-problems are solved by Gurobi as LPs.
The column-generation algorithm is implemented in the un-registered package [Linda.jl](https://github.com/mtanneau/Linda.jl).

We save the RMP every 10 iterations of the CG procedure, and at the final iteration.
This gives a dataset of structured LPs.

## Generating a sysimage

```bash
julia --trace-compile=exp/precompile.jl --project exp/snoop.jl
```

```julia
using PackageCompiler
PackageCompiler.create_sysimage([:JuMP, :LinearAlgebra, :SparseArrays, :Gurobi, :Linda, :SMPSReader]; project=".", sysimage_path="JuliaTSSP.so", precompile_statements_file="exp/precompile.jl");
```

## Generating RMP collection


1. Generate the list of jobs
    ```bash
    julia exp/jobs.jl > exp/jobs.txt
    ```

2. Run the jobs, e.g. with GNU `parallel`
    ```bash
    cat exp/jobs.txt | parallel -j1 --joblog colgen.log {}
    ```

3. The Master problems are located in `data/rmp`, in `.mps.bz2` format.
The naming convention is `<instance_name>_<scenarios>_<cg_iter>.mps.bz2`.