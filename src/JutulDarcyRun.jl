module JutulDarcyRun

using ArgParse
using DelimitedFiles: writedlm
using FileWatching.Pidfile: mkpidlock
using Jutul: Cells, JacobiPreconditioner, JutulCase, MultiModel, ParallelCSRContext,
    degrees_of_freedom_per_entity, get_1d_interpolator, get_variable, is_cell_major,
    jutul_message, matrix_layout, si_unit
using Jutul: export_mesh_vtu, physical_representation
using JutulDarcy: setup_case_from_afi, setup_case_from_data_file, simulate_reservoir,
    write_standard_output
import Jutul
import JutulDarcy

const MAX_PRECOMPILE_COMPONENTS = 5

include("cli.jl")
include("logging.jl")
include("simulation.jl")
include("restart_wells.jl")
include("precompile.jl")

end
