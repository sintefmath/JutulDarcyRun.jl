module JutulDarcyRunCUDA

using JutulDarcyRun: precompile_simulations, run_cli
import CUDA
using PrecompileTools: @compile_workload

device_available() = CUDA.functional() && CUDA.has_cuda_gpu()

function check_device()
    device_available() ||
        throw(ArgumentError("No functional CUDA device found. jdrun-cuda requires a CUDA-capable GPU and driver."))
    return nothing
end

function (@main)(args)
    return run_cli(args, "jdrun-cuda", :ka_cuda; check_device)
end

@compile_workload begin
    if device_available()
        precompile_simulations(:ka_cuda, "jdrun-cuda"; strict = false)
    end
end
precompile(main, (Vector{String},))

end
