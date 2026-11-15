# JutulDarcyRun.jl

Run JutulDarcy reservoir simulations from the terminal. Requires Julia 1.12 or later.

## Install

Start Julia with `julia -O3 --threads=auto`. Press `]` to open the package prompt.

```text
app add https://github.com/sintefmath/JutulDarcyRun.jl
```

For a GPU launcher, press Backspace at an empty package prompt and choose
the install command for your GPU.

```julia
using Pkg
Pkg.Apps.add(url="https://github.com/sintefmath/JutulDarcyRun.jl.git", subdir="apps/cuda")
Pkg.Apps.add(url="https://github.com/sintefmath/JutulDarcyRun.jl.git", subdir="apps/amdgpu")
```

Add the `bin` directory of your first Julia depot to `PATH`. The default is
`~/.julia/bin` on Linux and macOS and `%USERPROFILE%\.julia\bin` on Windows.

## Run a case

From the terminal, run

```sh
jdrun MODEL.DATA
jdrun MODEL.DATA --number-of-steps 1 --output-path results
jdrun --help
```

Inputs can be Eclipse `.DATA` decks, `.AFI` files, or Julia scripts that define a
top-level `case::JutulCase`. The CPU app uses KernelAbstractions by default.
Use `--backend cpu` or `--backend hypre` for the CPU path with HYPRE as the CPR
pressure solver. For GPU runs, use `jdrun-cuda` or `jdrun-amdgpu`
with a working CUDA or ROCm setup.

Results go in `<case>_jutul_output` next to the input. The Eclipse output includes
the grid (`.EGRID`), report-step states (`.UNRST`), and summary (`.SMSPEC` and
`.UNSMRY`). Add `--vtk` for ParaView files or `--csv` for summary time series.
`--output-substates` includes internal timesteps in the summary, CSV, and VTK
output for single-process runs.

### Solver options

`--precond cpr` and `--precond cprw` select a two-stage preconditioner. Use
`--global-smoother` to choose its full-system stage. `--precond ilu0`,
`--precond dilu`, and `--precond spai0` select standalone smoothers. Jacobi is
available on CPU, and `--global-smoother vendor_ilu` uses the GPU library's ILU.
MPI supports `cpr` for pressure preconditioning.

`--linear-solver` selects GMRES, FGMRES, or BiCGSTAB. The GPU launchers use
Float32 FGMRES by default. Select `--single-precision-linear-solver false` to
use Float64, or choose another Krylov solver. The `ka` CPU backend also accepts
Float32 FGMRES.

`--timesteps none` disables the adaptive selectors. Initial and maximum step
sizes are in days, and timestep cuts can still occur. The saturation-change
target is a fraction and applies to automatic timestepping with the HYPRE CPU
backend.

Pass Julia flags before `--` to override the launcher's thread count.

```sh
jdrun --threads=8 -- MODEL.DATA --backend ka --global-smoother dilu
```

### MPI

Use the `mpiexec` that matches your MPI.jl setup. All ranks need access to the input
and output paths. For Julia case scripts, set `split_wells = true` when building
the model. [MPI.jl configuration](https://juliaparallel.org/MPI.jl/stable/configuration/)

```sh
mpiexec -n 4 jdrun MODEL.DATA --backend mpi
```

On Windows, use `mpiexec -n 4 cmd.exe /c jdrun MODEL.DATA --backend mpi`.

## Develop locally

Clone the repository and start Julia from its root.

```sh
git clone https://github.com/sintefmath/JutulDarcyRun.jl
cd JutulDarcyRun.jl
julia -O3 --threads=auto
```

Use the package prompt (`]`) to run from your working copy.

```text
app develop .
```

### GPU launchers

Start Julia with the same flags from `apps/cuda` or `apps/amdgpu`. At the
package prompt, point the GPU project at the local runner and install the launcher.

```text
activate .
develop ../..
app develop .
```
