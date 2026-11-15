module JutulDarcyRunCUDAExt

import CUDA
using JutulDarcyRun: precompile_simulations
using PrecompileTools: @compile_workload

@compile_workload begin
    if CUDA.functional() && CUDA.has_cuda_gpu()
        precompile_simulations(:ka_cuda; strict = false)
    end
end

end
