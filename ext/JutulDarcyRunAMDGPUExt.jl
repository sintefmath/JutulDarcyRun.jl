module JutulDarcyRunAMDGPUExt

import AMDGPU
using JutulDarcyRun: precompile_simulations
using PrecompileTools: @compile_workload

@compile_workload begin
    if AMDGPU.functional() && AMDGPU.has_rocm_gpu()
        precompile_simulations(:ka_amd; strict = false)
    end
end

end
