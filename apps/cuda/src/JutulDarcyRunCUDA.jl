module JutulDarcyRunCUDA

using JutulDarcyRun: run_cli

function check_device()
    @eval import CUDA
    Base.invokelatest(() -> CUDA.functional() && CUDA.has_cuda_gpu()) ||
        throw(ArgumentError("No functional CUDA device found. jdrun-cuda requires a CUDA-capable GPU and driver."))
    return nothing
end

function (@main)(args)
    return run_cli(args, "jdrun-cuda", :ka_cuda; check_device)
end

end
