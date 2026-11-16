module JutulDarcyRunMPI

import HYPRE
import MPI
import PartitionedArrays
using Jutul: Cells, JutulCase, degrees_of_freedom_per_entity, simulate!
import JutulDarcy
using JutulDarcy: ReservoirSimResult, setup_mini_wellcase, setup_reservoir_simulator
using JutulDarcyRun: MAX_PRECOMPILE_COMPONENTS, check_solver_layout, limit_case, load_case,
    normalize_reservoir_variables!, parse_commandline, precompile_compositional_case,
    precompile_olympus_setup, precompile_smoothers, publish_output, run_cli, simulation_kwargs, with_run_log
using PrecompileTools: @compile_workload

include("mpi.jl")

function (@main)(args)
    return run_cli(args, "jdrun-mpi", :mpi; run_case = run_mpi_case)
end
julia_main()::Cint = main(ARGS)

function can_run_mpi_precompile()
    MPI.Finalized() && return false
    MPI.Initialized() && return MPI.Comm_size(MPI.COMM_WORLD) == 1
    # An uninitialized compiler worker can inherit a rank's launcher environment.
    return !any(key -> haskey(ENV, key), (
        "OMPI_COMM_WORLD_RANK", "PMI_RANK", "PMIX_RANK", "MPI_LOCALRANKID",
        "MV2_COMM_WORLD_RANK", "SLURM_PROCID",
    ))
end

function precompile_mpi()
    can_run_mpi_precompile() || return nothing
    MPI.Init()
    comm = MPI.COMM_WORLD
    MPI.Comm_size(comm) == 1 || return nothing
    precompile_olympus_setup(; split_wells = true)
    args = parse_commandline(["case.jl", "--info-level", "-1"], "jdrun-mpi", :mpi)
    mktempdir() do dir
        for physics in (:single_phase, :immiscible_2ph)
            case = setup_mini_wellcase(physics; dims = (4, 1, 1), nstep = 1,
                total_time = 60.0, block_backend = true, split_wells = true)
            run_mpi_simulation(case, dir, "case", comm, args)
        end
        for ncomponents in 2:MAX_PRECOMPILE_COMPONENTS
            case = precompile_compositional_case(ncomponents; split_wells = true)
            run_mpi_simulation(case, dir, "case", comm, args)
        end
        for ncomponents in 1:MAX_PRECOMPILE_COMPONENTS
            precompile_smoothers(Val(ncomponents), :default)
        end
    end
    return nothing
end

@compile_workload begin
    precompile_mpi()
end
precompile(main, (Vector{String},))

end
