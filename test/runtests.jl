using Jutul, JutulDarcy, JutulDarcyRun, Test

function options(mode, flags...)
    args = JutulDarcyRun.parse_commandline(["case.jl", "--info-level", "-1", flags...], "jdrun-ka", mode)
    return JutulDarcyRun.simulation_kwargs(args, mode)
end

@testset "Command-line help and verbosity" begin
    app = JutulDarcyRun
    settings = app.commandline_settings("jdrun-ka", :ka)
    help = sprint(io -> app.ArgParse.show_help(io, settings; exit_when_done = false))
    words = join(split(help), " ")
    @test startswith(words, "usage: jdrun-ka [options] filename")
    @test occursin("Default is 0.05.", words)
    @test occursin("Defaults to partial_sparsity with zero rebuild tolerance and memory otherwise.", words)
    @test occursin("type: Float64, default: 1.0", words)
    for mode in (:ka, :default, :mpi, :ka_cuda, :ka_amd)
        args = app.parse_commandline(["case.jl"], "jdrun-ka", mode)
        @test args[:info_level] === 1.0
        @test args[:amg_rebuild_tolerance] === nothing
        @test options(mode, "--info-level", "0.5")[:info_level] === 0.5
    end
end

@testset "Quiet deck setup during precompilation" begin
    app = JutulDarcyRun
    deck = joinpath(pkgdir(JutulDarcy), "test", "data", "thermal", "TRACER_SMALL.DATA")
    mktempdir() do dir
        output = joinpath(dir, "setup.log")
        redirect_stdio(; stdout = output, stderr = output) do
            app.load_case(deck, ".data"; verbose = false, split_wells = false)
        end
        @test occursin("RUNSPEC/TRACERS", read(output, String))
        case = redirect_stdio(; stdout = output, stderr = output) do
            app.load_precompile_deck(deck)
        end
        @test case isa JutulCase
        @test isempty(read(output, String))
    end
end

@testset "CPRW default and override" begin
    for mode in (:ka, :default, :ka_cuda, :ka_amd)
        @test options(mode)[:precond] == :cprw
        @test options(mode)[:linear_solver_arg].smoother_type == :ka_dilu
        @test options(mode, "--precond", "cpr")[:precond] == :cpr
    end
    @test !haskey(options(:mpi), :precond)
    @test options(:mpi, "--precond", "cpr")[:precond] == :cpr
    case = JutulDarcy.setup_mini_wellcase(:immiscible_2ph; nstep = 1)
    for mode in (:ka, :default)
        _, config = setup_reservoir_simulator(case; options(mode)...)
        @test config[:linear_solver].preconditioner.variant == :cprw
    end
    single = JutulDarcy.setup_mini_wellcase(:single_phase; nstep = 1, block_backend = true)
    args = JutulDarcyRun.parse_commandline(["case.jl"], "jdrun-ka", :ka)
    @test JutulDarcyRun.case_preconditioner(case, args, :cprw) == :cprw
    @test JutulDarcyRun.case_preconditioner(single, args, :cprw) == :dilu
    tracer = JutulDarcyRun.precompile_case(:deck_tracer, false)
    @test JutulDarcyRun.case_preconditioner(tracer, args, :cprw) == :cpr
    csc = JutulDarcy.setup_mini_wellcase(:immiscible_2ph; nstep = 1, backend = :csc)
    @test JutulDarcyRun.case_preconditioner(csc, args, :cprw) == :cpr
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
    for mode in (:ka_cuda, :ka_amd)
        @test options(mode)[:wells_on_device]
        @test !options(mode, "--wells-on-device", "false")[:wells_on_device]
    end
    @test !haskey(options(:ka), :wells_on_device)
    @test options(:ka, "--ka-workgroupsize", "128")[:ka_workgroupsize] == 128
end

@testset "CPR update settings" begin
    case = JutulDarcy.setup_mini_wellcase(:immiscible_2ph; nstep = 1)
    defaults = options(:ka)[:linear_solver_arg]
    @test keys(defaults) == (:smoother_type,)
    _, config = setup_reservoir_simulator(case; options(:ka, "--amg-rebuild-tolerance", "0")...)
    @test config[:linear_solver].preconditioner.pressure_precond.reuse == :partial_sparsity

    kwarg = options(:ka, "--amg-rebuild-tolerance", "0.1", "--update-type", "sparsity",
        "--update-type-partial", "partial_operators", "--update-interval", "step",
        "--update-interval-partial", "ministep", "--partial-update", "false")
    _, config = setup_reservoir_simulator(case; kwarg...)
    cpr = config[:linear_solver].preconditioner
    @test cpr.amg_rebuild_tolerance == 0.1
    @test cpr.pressure_precond.reuse == :sparsity
    @test cpr.pressure_precond.reuse_partial == :partial_operators
    @test cpr.update_interval == :step
    @test cpr.update_interval_partial == :ministep
    @test !cpr.partial_update
end

@testset "Case normalization" begin
    case = JutulDarcy.setup_mini_wellcase(:immiscible_2ph; nstep = 1)
    reservoir = JutulDarcy.reservoir_model(case)
    parameters = case.parameters[:Reservoir]
    parameters[:FluidVolume] .*= 1.25
    volume = copy(parameters[:FluidVolume])
    JutulDarcyRun.normalize_reservoir_variables!(case)
    @test parameters[:StaticFluidVolume] == volume
    @test case.state0[:Reservoir][:FluidVolume] == volume
    @test all(iszero, reservoir.secondary_variables[:FluidVolume].expansion)
    @test all(iszero, case.state0[:Reservoir][:CapillaryPressure])
    pc = reservoir.secondary_variables[:CapillaryPressure]
    @test only(only(pc.pc))(0.4) == 0.0
    JutulDarcyRun.normalize_reservoir_variables!(case)
    @test reservoir.secondary_variables[:CapillaryPressure] === pc

    olympus = JutulDarcyRun.precompile_case(:deck_olympus, false)
    reservoir = JutulDarcy.reservoir_model(olympus)
    rock = reservoir.secondary_variables[:FluidVolume]
    JutulDarcyRun.normalize_reservoir_variables!(olympus)
    @test reservoir.secondary_variables[:FluidVolume] === rock
    @test length(olympus.dt) == 1
end

@testset "Restart well completions" begin
    app = JutulDarcyRun
    geo = JutulDarcy.GeoEnergyIO
    mktempdir() do dir
        deck = read(joinpath(@__DIR__, "data", "OLYMPUS_SMALL.DATA"), String)
        path = joinpath(dir, "case.DATA")
        write(path, replace(deck, "\nPROPS" => "\nACTNUM\n 1 1 0 1 /\nPROPS"))
        case = app.load_case(path, ".data"; verbose = false, split_wells = false)
        (; dims, cell_map) = app.grid_layout(case)
        @test cell_map == [1, 2, 4]
        pressure(step) = Float32[10, 20, 30] .* step
        writer = geo.StandardIO.output
        blocks = map(1:2) do step
            headers = writer.standard_headers(dims, length(cell_map); elapsed = step * 86400, report_step = step)
            writer.restart_block(headers, (geo.BinaryKeyword("PRESSURE", pressure(step)),); report_step = step)
        end
        for model_case in (case, JutulDarcy.merge_similar_wells(case))
            restart = geo.write_restart(blocks, joinpath(dir, "CASE"))
            app.add_restart_wells(restart, model_case)
            app.add_restart_wells(restart, model_case)
            reports, raw = geo.read_restart(restart; extra_out = true)
            @test length(reports) == 2
            for (step, report) in enumerate(reports)
                @test report["PRESSURE"] == pressure(step)
                @test report["days"] == step
                @test length(geo.keyword_values(raw[step], "ZWEL")) == 1
                nwells, ncwmax, niwelz, nzwelz, niconz = report["INTEHEAD"][[17, 18, 25, 28, 33]]
                names = reshape(report["ZWEL"], nzwelz, nwells)
                counts = reshape(report["IWEL"], niwelz, nwells)
                connections = reshape(report["ICON"], niconz, ncwmax, nwells)
                @test counts[5, :] == [1, 1]
                @test Dict(names[1, i] => Tuple(connections[2:4, 1, i]) for i in 1:nwells) ==
                    Dict("I1" => (1, 1, 1), "P1" => (4, 1, 1))
            end
        end
    end
end

@testset "Run log" begin
    mktempdir() do dir
        terminal = joinpath(dir, "terminal.txt")
        redirect_stdio(; stdout = terminal, stderr = terminal) do
            JutulDarcyRun.with_run_log(dir, "case", "input.DATA") do
                println("Report step 1")
                println(stderr, "Case diagnostic")
            end
            @test_throws ArgumentError JutulDarcyRun.with_run_log(dir, "failed", "input.DATA") do
                throw(ArgumentError("Invalid input"))
            end
        end
        log = read(joinpath(dir, "case.DBG"), String)
        output = read(terminal, String)
        for message in ("Report step 1", "Case diagnostic")
            @test occursin(message, log)
            @test occursin(message, output)
        end
        @test occursin("Run complete", log)
        @test occursin("Invalid input", read(joinpath(dir, "failed.DBG"), String))
    end
end

@testset "Krylov layout checks use the reservoir" begin
    args = JutulDarcyRun.parse_commandline(["case.jl", "--linear-solver", "gmres"], "jdrun-mpi", :mpi)
    block = JutulDarcy.setup_mini_wellcase(:immiscible_2ph; nstep = 1, split_wells = true)
    scalar = JutulDarcy.setup_mini_wellcase(:immiscible_2ph; nstep = 1, block_backend = false)
    @test JutulDarcyRun.check_solver_layout(block, :mpi, args) === nothing
    @test_throws ArgumentError JutulDarcyRun.check_solver_layout(scalar, :default, args)
end
