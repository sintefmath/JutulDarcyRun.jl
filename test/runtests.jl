using Jutul, JutulDarcy, JutulDarcyRun, Test

function options(mode, flags...)
    launcher = mode in (:ka_cuda, :ka_amd) ? mode : :cpu
    args = JutulDarcyRun.parse_commandline(["case.jl", "--info-level", "-1", flags...], "jdrun", launcher)
    return JutulDarcyRun.simulation_kwargs(args, mode)
end

@testset "Standalone preconditioners" begin
    case = JutulDarcy.setup_mini_wellcase(:immiscible_2ph; nstep = 1)
    for mode in (:ka, :default)
        for (name, expected) in (
                ("ilu0", Jutul.ILUZeroPreconditioner),
                ("jacobi", Jutul.JacobiPreconditioner),
                ("dilu", Jutul.KASmootherPreconditioner),
                ("spai0", mode == :ka ? Jutul.KASmootherPreconditioner : Jutul.SPAI0Preconditioner),
            )
            _, config = setup_reservoir_simulator(case; options(mode, "--precond", name)...)
            preconditioner = config[:linear_solver].preconditioner
            @test preconditioner isa expected
            if name == "dilu"
                @test preconditioner.config isa Jutul.KAPreconditioners.DILU
            elseif name == "spai0" && mode == :ka
                @test preconditioner.config isa Jutul.KAPreconditioners.SPAI0
            end
        end
    end
end

@testset "Backend-specific options" begin
    @test_throws ArgumentError options(:mpi, "--precond", "cprw")
    @test_throws ArgumentError options(:mpi, "--output-substates")
    @test_throws ArgumentError options(:ka, "--global-smoother", "vendor_ilu")
    @test_throws ArgumentError options(:ka, "--precond", "jacobi", "--global-smoother", "dilu")
    @test_throws ArgumentError options(:default, "--timesteps", "none", "--timesteps-target-iterations", "8")
    @test_throws ArgumentError options(:default, "--timesteps", "iteration", "--timesteps-target-ds", "0.05")
    @test_throws ArgumentError options(:ka, "--linear-solver", "gmres", "--single-precision-linear-solver", "true")
    @test options(:default, "--timesteps-target-ds", "0.05")[:target_ds] == 0.05
    @test options(:ka_cuda)[:single_precision_linear_solver]
    @test options(:ka_cuda, "--single-precision-linear-solver", "false")[:linear_solver] == :fgmres
    @test !haskey(options(:ka_cuda, "--linear-solver", "gmres"), :single_precision_linear_solver)
end

@testset "Krylov layout checks use the reservoir" begin
    args = JutulDarcyRun.parse_commandline(["case.jl", "--linear-solver", "gmres"], "jdrun", :cpu)
    block = JutulDarcy.setup_mini_wellcase(:immiscible_2ph; nstep = 1, split_wells = true)
    scalar = JutulDarcy.setup_mini_wellcase(:immiscible_2ph; nstep = 1, block_backend = false)
    @test JutulDarcyRun.check_solver_layout(block, :mpi, args) === nothing
    @test_throws ArgumentError JutulDarcyRun.check_solver_layout(scalar, :default, args)
end
