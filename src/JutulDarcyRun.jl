module JutulDarcyRun

using ArgParse
using FileWatching.Pidfile: mkpidlock
using Jutul: JutulCase, jutul_message, si_unit
using JutulDarcy: setup_case_from_afi, setup_case_from_data_file, simulate_reservoir,
    write_standard_output
using PrecompileTools: @compile_workload, @setup_workload
import JutulDarcy

include("cli.jl")
include("simulation.jl")

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

@setup_workload begin
    case = JutulDarcy.setup_mini_wellcase(:immiscible_2ph; nstep = 1, total_time = si_unit(:day))
    @compile_workload begin
        mktempdir() do dir
            args = parse_commandline(["case.jl", "--info-level", "-1"], "jdrun", :cpu)
            run_simulation(case, dir, "case", :ka, args)
        end
    end
end

end
