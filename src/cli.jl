function parse_commandline(argv, command, backend)
    settings = ArgParseSettings(
        prog = command,
        description = "Run a JutulDarcy reservoir simulation and export unified Eclipse results.",
        version = string(pkgversion(@__MODULE__)),
        add_version = true,
        autofix_names = true,
    )
    @add_arg_table! settings begin
        "filename"
            help = "Eclipse .DATA deck, .AFI file, or Julia .jl script defining a JutulCase named case"
            required = true
    end

    add_arg_group!(settings, "Linear and nonlinear solvers")
    if backend == :cpu
        @add_arg_table! settings begin
            "--backend"
                help = "CPU solver backend (ka, hypre or mpi)"
                arg_type = Symbol
                default = :ka
                range_tester = x -> x in (:ka, :hypre, :mpi)
            "--method"
                help = "Nonlinear method (newton or nldd, which requires hypre or mpi)"
                arg_type = Symbol
                default = :newton
                range_tester = x -> x in (:newton, :nldd)
        end
    end
    @add_arg_table! settings begin
        "--linear-solver"
            help = "Krylov solver (gmres, fgmres or bicgstab)"
            arg_type = Symbol
            range_tester = x -> x in (:gmres, :fgmres, :bicgstab)
        "--precond"
            help = "Preconditioner (cprw, cpr, ilu0, jacobi or spai0)"
            arg_type = Symbol
            range_tester = x -> x in (:cprw, :cpr, :ilu0, :jacobi, :spai0)
        "--global-smoother", "--smoother-type"
            help = "Full-system smoother, such as ka_dilu, ka_ilu0, ka_spai0 or vendor_ilu (GPU)"
            arg_type = Symbol
        "--rtol"
            help = "Relative linear solver tolerance"
            arg_type = Float64
        "--tol-cnv"
            help = "CNV tolerance for mass balance"
            arg_type = Float64
        "--tol-mb"
            help = "Scaled L1 mass balance tolerance"
            arg_type = Float64
        "--tol-dp-well"
            help = "Well connection pressure drop tolerance in bar"
            arg_type = Float64
        "--tol-factor-final-iteration"
            help = "Factor applied to tolerances on the final Newton iteration"
            arg_type = Float64
        "--relaxation"
            help = "Use nonlinear solver damping (true or false)"
            arg_type = Bool
        "--max-nonlinear-iterations"
            help = "Maximum Newton iterations before cutting a timestep"
            arg_type = Int
            range_tester = x -> x > 0
        "--max-timestep-cuts"
            help = "Maximum timestep cuts per report step"
            arg_type = Int
            range_tester = x -> x >= 0
        "--presolve-wells"
            help = "Solve wells before each linearization (true or false)"
            arg_type = Bool
    end

    add_arg_group!(settings, "Timestepping")
    @add_arg_table! settings begin
        "--number-of-steps"
            help = "Limit the number of report steps"
            arg_type = Int
            range_tester = x -> x > 0
        "--timesteps"
            help = "Timestep selection (auto, iteration or none)"
            arg_type = Symbol
            range_tester = x -> x in (:auto, :iteration, :none)
        "--timesteps-target-iterations"
            help = "Target Newton iterations for timestep selection"
            arg_type = Int
            range_tester = x -> x > 0
        "--timesteps-target-ds"
            help = "Target saturation change for timestep selection (requires hypre)"
            arg_type = Float64
        "--initial-timestep"
            help = "Initial timestep in days"
            arg_type = Float64
            range_tester = x -> x > 0
        "--max-timestep"
            help = "Maximum timestep in days"
            arg_type = Float64
            range_tester = x -> x > 0
    end

    add_arg_group!(settings, "Input setup")
    @add_arg_table! settings begin
        "--skip-tracers"
            help = "Skip tracer transport in Eclipse .DATA decks"
            action = :store_true
    end

    add_arg_group!(settings, "Output and printing")
    @add_arg_table! settings begin
        "--output-path"
            help = "Output directory (defaults to <case>_jutul_output next to the input file)"
        "--vtk"
            help = "Also export reservoir states as a ParaView VTK time series"
            action = :store_true
        "--csv"
            help = "Also export field, well and group summary time series as CSV in SI units"
            action = :store_true
        "--output-substates"
            help = "Include internal timesteps in the Jutul results and summary"
            action = :store_true
        "--info-level"
            help = "Simulation verbosity (-1 for quiet)"
            arg_type = Int
            default = 1
        "--report-level"
            help = "Detail level of Jutul convergence reports"
            arg_type = Int
            default = 0
            range_tester = x -> x >= 0
        "--verbose"
            help = "Show additional input setup information"
            action = :store_true
        "--print-wells"
            help = "Print well results after simulation"
            action = :store_true
    end
    return parse_args(argv, settings; as_symbols = true)
end

function simulation_kwargs(args, mode)
    kwarg = Dict{Symbol, Any}(
        :mode => mode,
        :method => get(args, :method, :newton),
        :info_level => args[:info_level],
        :report_level => args[:report_level],
        :output_substates => args[:output_substates],
        :error_on_incomplete => true,
    )
    for key in (:linear_solver, :precond, :rtol, :tol_cnv, :tol_mb, :tol_dp_well,
            :tol_factor_final_iteration, :relaxation, :max_nonlinear_iterations,
            :max_timestep_cuts, :presolve_wells, :timesteps)
        value = args[key]
        isnothing(value) || (kwarg[key] = value)
    end
    for (key, value, scale) in (
            (:target_its, args[:timesteps_target_iterations], 1),
            (:target_ds, args[:timesteps_target_ds], 1),
            (:initial_dt, args[:initial_timestep], si_unit(:day)),
            (:max_dt, args[:max_timestep], si_unit(:day)))
        isnothing(value) || (kwarg[key] = value * scale)
    end
    linear_solver_arg = (;)
    if mode in (:default, :mpi)
        kwarg[:linear_solver_backend] = :cpu
        linear_solver_arg = (amg_type = :hypre,)
    end
    smoother_type = args[:global_smoother]
    if !isnothing(smoother_type)
        linear_solver_arg = (; linear_solver_arg..., smoother_type)
    end
    isempty(linear_solver_arg) || (kwarg[:linear_solver_arg] = linear_solver_arg)
    return kwarg
end

function run_cli(argv, command, backend; check_device = nothing)
    args = parse_commandline(argv, command, backend)
    try
        path = abspath(args[:filename])
        isfile(path) || throw(ArgumentError("Input file does not exist: $path"))
        name, extension = splitext(basename(path))
        extension = lowercase(extension)
        extension in (".data", ".afi", ".jl") ||
            throw(ArgumentError("Expected a .DATA, .AFI, or .jl input file, got $path"))

        mode = backend == :cpu ? args[:backend] : backend
        mode = mode == :hypre ? :default : mode
        method = get(args, :method, :newton)
        method == :newton || mode in (:default, :mpi) ||
            throw(ArgumentError("--method nldd requires --backend hypre or --backend mpi"))
        isnothing(args[:timesteps_target_ds]) || mode == :default ||
            throw(ArgumentError("--timesteps-target-ds requires --backend hypre"))
        !args[:skip_tracers] || extension == ".data" ||
            throw(ArgumentError("--skip-tracers applies to Eclipse .DATA input"))

        output_dir = args[:output_path]
        output_dir = isnothing(output_dir) ? joinpath(dirname(path), name * "_jutul_output") : abspath(output_dir)
        if mode == :mpi
            args[:output_substates] && throw(ArgumentError("MPI output supports report steps only"))
            @eval import PartitionedArrays
            Base.invokelatest(run_mpi_case, path, extension, output_dir, name, args)
            return 0
        end

        isnothing(check_device) || check_device()
        case = load_case(path, extension; verbose = args[:verbose], split_wells = method == :nldd,
            include_tracers = !args[:skip_tracers])
        # GPU extensions and Julia scripts can define methods after this call began.
        Base.invokelatest(run_simulation, case, output_dir, name, mode, args)
    catch err
        if err isa ArgumentError && !args[:verbose]
            println(stderr, command, ": ", err.msg)
            return 1
        end
        rethrow()
    end
    return 0
end
