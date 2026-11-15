module JutulDarcyRunAMDGPU

using JutulDarcyRun: precompile_simulations, run_cli
import AMDGPU
using PrecompileTools: @compile_workload

device_available() = AMDGPU.functional() && AMDGPU.has_rocm_gpu()

function check_device()
    device_available() ||
        throw(ArgumentError("No functional ROCm device found. jdrun-amdgpu requires a ROCm-compatible GPU and driver."))
    return nothing
end

function (@main)(args)
    return run_cli(args, "jdrun-amdgpu", :ka_amd; check_device)
end

@compile_workload begin
    if device_available()
        precompile_simulations(:ka_amd, "jdrun-amdgpu"; strict = false)
    end
end
precompile(main, (Vector{String},))

end
