module JutulDarcyRunAMDGPU

using JutulDarcyRun: run_cli

function check_device()
    @eval import AMDGPU
    Base.invokelatest(() -> AMDGPU.functional() && AMDGPU.has_rocm_gpu()) ||
        throw(ArgumentError("No functional ROCm device found. jdrun-amdgpu requires a ROCm-compatible GPU and driver."))
    return nothing
end

function (@main)(args)
    return run_cli(args, "jdrun-amdgpu", :ka_amd; check_device)
end

end
