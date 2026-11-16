module JutulDarcyRunHYPRE

import HYPRE
using JutulDarcyRun: precompile_simulations, run_cli
using PrecompileTools: @compile_workload

function (@main)(args)
    return run_cli(args, "jdrun-hypre", :default)
end
julia_main()::Cint = main(ARGS)

@compile_workload begin
    precompile_simulations(:default, "jdrun-hypre")
end
precompile(main, (Vector{String},))

end
