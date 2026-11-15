# JutulDarcyRun.jl

Run JutulDarcy reservoir simulations from the terminal with Julia 1.12 or later.

## Install

From a local checkout, run this at Julia's package prompt (`]`).

```text
app develop /path/to/JutulDarcyRun.jl
```

Use `app add /path/to/JutulDarcyRun.jl` to install from a committed checkout.
Add the `bin` subdirectory of your first Julia depot to `PATH`. With the default
depot, this is `~/.julia/bin` on Linux and macOS or `%USERPROFILE%\.julia\bin` on Windows.
See the [Pkg apps documentation](https://pkgdocs.julialang.org/v1/apps/) for more installation options.

## Run

```text
jdrun MODEL.DATA
jdrun MODEL.DATA --backend hypre --number-of-steps 10
jdrun-cuda MODEL.AFI --output-path results
jdrun-amdgpu case.jl
jdrun --threads=8 -- MODEL.DATA
```

`jdrun` uses KernelAbstractions on the CPU by default. `--backend hypre` selects
the CPU solver with HYPRE for the CPR pressure solve. The GPU commands require a
functional CUDA or ROCm device at runtime. All three commands can be installed
on a CPU-only machine and start Julia with `-O3 --threads=auto` by default.

Each command accepts Eclipse `.DATA` decks, `.AFI` files with RESQML grids, and
Julia `.jl` scripts. A script runs in Julia's main workspace and
must define a top-level `case::JutulCase`. This example builds a case from a
deck next to the script.

```julia
using JutulDarcy

case = setup_case_from_data_file(joinpath(@__DIR__, "MODEL.DATA"))
```

Use `--help` on any command for its solver, tolerance, timestepping, and output
options. `--initial-timestep` and `--max-timestep` take values in days. Omitted
solver options use JutulDarcy's defaults, including Eclipse `TUNING` where supported.

## Output

The default output directory is `<case>_jutul_output` next to the input file.
Jutul writes native report steps to temporary files during the run. After a
successful solve, JutulDarcy's StandardIO integration consolidates them into
`<case>.SMSPEC`, `<case>.UNSMRY`, `<case>.UNRST`, and `<case>.EGRID`.
Subsequent runs replace the standard files for the same case. The writer uses
deck grid data when present and a synthetic grid for cases without it.
