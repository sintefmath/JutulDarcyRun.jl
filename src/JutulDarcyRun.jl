module JutulDarcyRun

using ArgParse
using DelimitedFiles: writedlm
using FileWatching.Pidfile: mkpidlock
using Jutul: Cells, JacobiPreconditioner, JutulCase, MultiModel, ParallelCSRContext,
    degrees_of_freedom_per_entity, is_cell_major, jutul_message, matrix_layout, si_unit, simulate!
using Jutul: export_mesh_vtu, physical_representation
using JutulDarcy: setup_case_from_afi, setup_case_from_data_file, simulate_reservoir,
    write_standard_output, setup_reservoir_simulator, ReservoirSimResult
using PrecompileTools: @compile_workload
import JutulDarcy
import MPI

include("cli.jl")
include("simulation.jl")
include("mpi.jl")
include("precompile.jl")

function (@main)(args)
    return run_cli(args, "jdrun", :cpu)
end

@compile_workload begin
    precompile_simulations(:ka)
    # HYPRE also initializes MPI when used as a single-process preconditioner.
    if can_run_mpi_precompile()
        precompile_simulations(:default)
    end
    precompile(run_mpi_case, (String, String, String, String, Dict{Symbol, Any}))
end

end
