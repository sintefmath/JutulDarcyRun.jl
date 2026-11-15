module JutulDarcyRunMPIExt

import MPI
import PartitionedArrays
using JutulDarcy: setup_mini_wellcase
using JutulDarcyRun: can_run_mpi_precompile, parse_commandline, run_mpi_simulation
using PrecompileTools: @compile_workload

@compile_workload begin
    if can_run_mpi_precompile()
        MPI.Init()
        comm = MPI.COMM_WORLD
        if MPI.Comm_size(comm) == 1
            case = setup_mini_wellcase(:immiscible_2ph;
                dims = (4, 1, 1), nstep = 1, total_time = 60.0,
                block_backend = true, split_wells = true)
            args = parse_commandline(["case.jl", "--backend", "mpi", "--info-level", "-1"], "jdrun", :cpu)
            mktempdir() do dir
                run_mpi_simulation(case, dir, "case", comm, args)
            end
        end
    end
end

end
