module JutulDarcyRunKA

using JutulDarcyRun: precompile_simulations, run_cli
using PrecompileTools: @compile_workload

function (@main)(args)
    return run_cli(args, "jdrun-ka", :ka)
end

@compile_workload begin
    precompile_simulations(:ka, "jdrun-ka")
end
precompile(main, (Vector{String},))

end
