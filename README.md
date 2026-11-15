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
jdrun MODEL.DATA --global-smoother ka_dilu
jdrun MODEL.DATA --skip-tracers
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
`--skip-tracers` disables tracer transport when setting up an Eclipse `.DATA` deck.
`--global-smoother` (alias `--smoother-type`) selects the full-system smoother.
For example, `ka_dilu`, `ka_ilu0`, and `ka_spai0` select KA implementations,
and `vendor_ilu` selects the native GPU implementation on CUDA and AMDGPU.
`--precond` still selects the preconditioner, such as CPR.

## MPI

Use the `mpiexec` for the MPI implementation configured in the app's Julia
environment. Input and output paths must be accessible to all ranks. See the
[MPI.jl configuration guide](https://juliaparallel.org/MPI.jl/stable/configuration/).

```sh
mpiexec -n 4 jdrun MODEL.DATA --backend mpi --output-path results
```

On Windows, launch the Pkg-generated batch file through `cmd.exe`.

```powershell
mpiexec -n 4 cmd.exe /c jdrun MODEL.DATA --backend mpi --output-path results
```

Each rank loads the input and runs its part of the solve. For `.jl` inputs,
construct `case` with `split_wells = true` in the model setup on every rank.
Rank zero writes the unified Eclipse files after Jutul consolidates the
distributed report steps. MPI defaults to CPR for multiphase mass-only
systems and to ILU(0) for scalar systems or systems with extra equations,
such as energy or tracer transport.

## Precompilation

The app precompiles one-step cases for CPU KA and HYPRE, covering single-phase,
immiscible, black-oil, compositional, geothermal, CO₂/brine, gas/water, and
Eclipse tracer transport. The workload includes thermal variants and exercises
native Jutul output followed by the standard export. CUDA and AMDGPU extensions
run the same workload when their device is functional during precompilation.
The MPI workload runs a four-cell, one-step CPR/HYPRE solve through the
partitioned simulator on a single rank outside an MPI job.

## Output

The default output directory is `<case>_jutul_output` next to the input file.
Jutul writes native report steps to temporary files during the run. After a
successful solve, JutulDarcy's StandardIO integration consolidates them into
`<case>.SMSPEC`, `<case>.UNSMRY`, `<case>.UNRST`, and `<case>.EGRID`.
Subsequent runs replace the standard files for the same case. The writer uses
deck grid data when present and a synthetic grid for cases without it.

```text
jdrun MODEL.DATA --vtk --csv
```

`--vtk` writes `<case>.pvd` and a `.vtu` file per saved reservoir state using
the simulation mesh. Open the `.pvd` file in ParaView to browse the states by
step number. `--csv` writes `<case>.csv` with field, well, and group summary
values in SI units. Each row contains the elapsed time in seconds, scope, name,
summary keyword, and value. With MPI, rank zero writes both optional exports
after consolidation.
