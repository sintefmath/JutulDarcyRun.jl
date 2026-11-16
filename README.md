# JutulDarcyRun.jl

Run JutulDarcy reservoir simulations from the terminal. Requires Julia 1.12 or later.

## Install

Start Julia with `julia -O3 --threads=8` and run the line for the backend you
want.

```julia
using Pkg
repo = "https://github.com/sintefmath/JutulDarcyRun.jl.git"
Pkg.Apps.add(url=repo, subdir="apps/ka")      # jdrun-ka (CPU KernelAbstractions)
Pkg.Apps.add(url=repo, subdir="apps/hypre")   # jdrun-hypre
Pkg.Apps.add(url=repo, subdir="apps/mpi")     # jdrun-mpi
Pkg.Apps.add(url=repo, subdir="apps/cuda")    # jdrun-cuda
Pkg.Apps.add(url=repo, subdir="apps/amdgpu")  # jdrun-amdgpu
```

Add the `bin` directory of your first Julia depot to `PATH`. The default is
`~/.julia/bin` on Linux and macOS and `%USERPROFILE%\.julia\bin` on Windows.

## Run a case

From the terminal, run

```sh
jdrun-ka MODEL.DATA
jdrun-ka MODEL.DATA --number-of-steps 1 --output-path results
jdrun-ka --help
```

Inputs can be Eclipse `.DATA` decks, `.AFI` files, or Julia scripts that define a
top-level `case::JutulCase`.

Results go in `<case>_jutul_output` next to the input. The Eclipse output includes
the grid (`.EGRID`), report-step states and well completions for plotting
(`.UNRST`), and summary (`.SMSPEC` and `.UNSMRY`). The live `<case>.DBG` log
captures terminal output. For MPI, it records rank 0. Add `--vtk` for ParaView
files or `--csv` for summary time series.
`--output-substates` includes internal timesteps in the summary, CSV, and VTK
output for single-process runs.

### Solver options

DILU is the default smoother. Use `--global-smoother ilu0` or
`--global-smoother spai0` to change it in CPR or CPRW. `--precond dilu`,
`--precond ilu0`, and `--precond spai0` select standalone smoothers.

Use `--amg-rebuild-tolerance`, `--update-type`, and `--update-type-partial`
to control KA AMG rebuilds. `--update-interval` and `--update-interval-partial`
set the CPR update schedules. Omitted options use JutulDarcy's defaults.

`--linear-solver` selects GMRES, FGMRES, or BiCGSTAB. The GPU launchers use
Float32 FGMRES by default. Select `--single-precision-linear-solver false` to
use Float64, or choose another Krylov solver. The `ka` CPU backend also accepts
Float32 FGMRES. Use `--ka-workgroupsize` to set the workgroup size.

GPU launchers merge similar wells and evaluate them on the GPU by default. Set
`--wells-on-device false` to disable merging and evaluate wells on the CPU.
Cases with wells that cannot be merged, such as tracer decks, use CPU well evaluation.

`--timesteps none` disables the adaptive selectors. Initial and maximum step
sizes are in days, and timestep cuts can still occur. The saturation-change
target is a fraction and applies to automatic timestepping with the HYPRE CPU
backend.

### Thread count

Each launcher starts with eight Julia threads.
To change the count for one run, pass `--threads` before the first `--`.
Pass the case path and solver options after it. [Julia Pkg app flags](https://pkgdocs.julialang.org/v1/apps/#Runtime-Julia-Flags)

```sh
jdrun-ka --threads=4 -- MODEL.DATA --precond cpr
jdrun-cuda --threads=16 -- MODEL.DATA
```

For a different default in a local checkout, edit `julia_flags` in the
launcher's `Project.toml` and rerun `app develop .` from that project to
regenerate the launcher.

### MPI

Use the `mpiexec` that matches your MPI.jl setup. All ranks need access to the input
and output paths. For Julia case scripts, set `split_wells = true` when building
the model. [MPI.jl configuration](https://juliaparallel.org/MPI.jl/stable/configuration/)

```sh
mpiexec -n 4 jdrun-mpi MODEL.DATA
```

On Windows, use `mpiexec -n 4 cmd.exe /c jdrun-mpi MODEL.DATA`.
Pass `--threads=2 --` before the case path to set the thread count per rank.

## Develop locally

Clone the repository and start Julia from the variant you want to develop.

```sh
git clone https://github.com/sintefmath/JutulDarcyRun.jl
cd JutulDarcyRun.jl/apps/ka
julia -O3 --threads=8
```

At the package prompt (`]`), point the app at the shared working copy and
install the launcher.

```text
activate .
develop ../..
app develop .
```
