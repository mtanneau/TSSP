# TSSP

Experiments for Two-Stage Stochastic Programming instances

## Problem formulation

A TSSP problem writes
\begin{align}
    \min_{x, y} \ \ \ & c^{T}x + \sum_{i} p_{i} q_{i}^{T}y_{i}\\
    s.t. \ \ \
    & l^{(1)} \leq A x \leq u^{(1)}\\
    & l^{(2)}_{i} \leq T_{i} x + W_{i} y_{i} \leq u^{(2)}_{i} & \forall i\\
    & x, y_{i} \geq 0
\end{align}

We transform it into standard form (by adding appropriate slacks)
\begin{align}
    \min_{x, y} \ \ \ & c^{T}x + \sum_{i} p_{i} q_{i}^{T}y_{i}\\
    s.t. \ \ \
    & A x = b\\
    & T_{i} x + W_{i} y_{i} = h_{i}& \forall i\\
    & x, y_{i} \geq 0
\end{align}
whose dual writes
\begin{align}
    \min_{\eta, \theta} \ \ \ & -b^{T}\eta - \sum_{i} h_{i}^{T}\theta_{i}\\
    s.t. \ \ \
    & A^{T} \eta + \xi + \sum_{i} T_{i}^{T} \theta_{i} = c \\
    & W_{i}^{T} \theta_{i} \leq p_{i} q_{i}& \forall i\\
    & \xi \geq 0
\end{align}


## DW decomposition

Create one sub-problem per scenario

### Master problem

The Master problem writes

$$
    \begin{array}{rl}
    \displaystyle \min_{\eta, \lambda} \ \ \ & -b^{T} \eta + \sum_{i, \rho} \lambda_{i, \rho} c_{i, \rho} + \sum_{i, \omega} \lambda_{i, \omega} c_{i, \omega}\\
    s.t. \ \ \ 
    & A^{T} \eta + \xi + \sum_{i, \rho} \lambda_{i, \rho} t_{i, \rho} + \sum_{i, \omega} \lambda_{i, \omega} t_{i, \omega} = c\\
    & \sum_{\omega} \lambda_{i, \omega} = 1 & \forall i\\
    & \xi, \lambda \geq 0
    \end{array}
$$
where $c_{i, \omega} = -h_{i}^{T}\omega$ and $t_{i, \omega} = T_{i}^{T} \omega$.

Given that $\eta$ is a free variable, we will write it as the difference of two non-negative variables.

### Sub-problems
The $i$-th sub-problem writes

$$
    \min_{\theta} \ \ \ (h_{i} - T_{i} y)^{T} \theta\\
    s.t. \ \ \ W_{i}^{T} \theta \leq p_{i} q_{i}
$$

## Installation

1. Clone/download this repository
1. Download and install commercial solvers
    * CPLEX
    * Gurobi
    * Mosek
1. Install unregistered packages:
    * `SMPSReader`
    * `Linda` (on the `callback` branch)
    * `UnitBlockAngular`
1. Install registered packages
    ```bash
    julia --project -e 'using Pkg; Pkg.instantiate()'
    ```

## Instances

We consider the TSSP instances with at least 1000 scenarios.
We remove `assets-large` which is solved in a single CG iteration.

Each instance is solved by column-generation, with Gurobi (default parameters) as the MP solver.
Sub-problems are solved by Gurobi as LPs.

We save the RMP every 10 iterations of the CG procedure, and at the final iteration.
This gives a dataset of structured LPs.

| `.cor` file | `.tim` file | `.sto` file |
|:------------|:------------|:------------|
`env.cor`, `env.cor.diss` | `env.tim` | `env.sto.N`
`4node.cor`, `4node.cor.base` | `4node.tim` | `4node.sto.N`
`phone.cor` | `phone.tim` | `phone.sto`
`stormG2.cor` | `stormG2.tim` | `stormG2_1000.sto`

We remove `assets` instances because they are solved in a single iteration.

## Generating a sysimage

```bash
julia --trace-compile=exp/precompile.jl --project exp/snoop.jl
```

```julia
using PackageCompiler
PackageCompiler.create_sysimage([:JuMP, :LinearAlgebra, :SparseArrays, :Gurobi, :CPLEX, :Mosek, :MosekTools, :Linda, :SMPSReader]; project=".", sysimage_path="JuliaTSSP.so", precompile_statements_file="exp/precompile.jl");
```

## Running the experiments

Running the experiments is done in two steps.
First, one needs to run the column-generation algorithm to generate
the RMP instances.
Then, the saved MPS files are solved with various solvers.

### Column-generation

1. Generate the list of jobs
    ```bash
    julia exp/colgen/jobs.jl > exp/colgen/jobs.txt
    ```

2. Run the jobs, e.g. with GNU `parallel`
    ```bash
    cat exp/colgen/jobs.txt | parallel -j1 --joblog colgen.log {}
    ```

### IPM


