module JutulDarcyRun

using ArgParse
using FileWatching.Pidfile: mkpidlock
using Jutul: Cells, JutulCase, degrees_of_freedom_per_entity, jutul_message, si_unit, simulate!
using Jutul: get_simulator_model
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

module CUDAApp
using ..JutulDarcyRun: run_cli
function (@main)(args)
    return run_cli(args, "jdrun-cuda", :ka_cuda)
end
end

module AMDGPUApp
using ..JutulDarcyRun: run_cli
function (@main)(args)
    return run_cli(args, "jdrun-amdgpu", :ka_amd)
end
end

@compile_workload begin
    precompile_simulations(:ka)
    precompile_simulations(:default)
    # MPI initialization and collectives must run in the launched MPI job.
    precompile(run_mpi_case, (String, String, String, String, Dict{Symbol, Any}))
end

end
