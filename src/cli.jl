const STANDALONE_PRECONDITIONERS = (:ilu0, :dilu, :jacobi, :spai0)
const CPR_REBUILD_OPTIONS = (:amg_rebuild_tolerance, :update_type, :update_type_partial,
    :update_interval, :update_interval_partial, :partial_update)

function commandline_settings(command, mode)
    gpu = mode in (:ka_cuda, :ka_amd)
    choices(values) = "{" * join(values, ",") * "}"
    methods = (:newton, :nldd)
    linear_solvers = (:gmres, :fgmres, :bicgstab)
    preconditioners = (:cprw, :cpr, STANDALONE_PRECONDITIONERS...)
    global_smoothers = (:dilu, :ilu0, :ilu, :jacobi, :spai0, :vendor_ilu, :ka_ilu0, :ka_dilu, :ka_spai0)
    amg_update_types = (:memory, :sparsity, :partial_sparsity, :partial_operators, :operators)
    cpr_update_intervals = (:once, :iteration, :ministep, :step)
    timestep_selectors = (:auto, :iteration, :none)
    # Leave delegated options unset so case and solver defaults still apply.
    input_default(key) = JutulDarcy.reservoir_simulation_input_defaults(nothing, key)
    settings = ArgParseSettings(
        prog = command,
        usage = "usage: $command [options] filename",
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
    if mode in (:default, :mpi)
        @add_arg_table! settings begin
            "--method"
                help = "Nonlinear method"
                metavar = choices(methods)
                arg_type = Symbol
                default = :newton
                range_tester = x -> x in methods
        end
    end
    @add_arg_table! settings begin
        "--linear-solver"
            help = gpu ? "Krylov solver" : "Krylov solver. Default is gmres."
            metavar = choices(linear_solvers)
            arg_type = Symbol
            default = gpu ? :fgmres : nothing
            range_tester = x -> x in linear_solvers
        "--single-precision-linear-solver"
            help = gpu ? "Use Float32. Defaults to true with FGMRES and false with other solvers." :
                "Use Float32 with FGMRES on CPU KA"
            metavar = "{true,false}"
            arg_type = Bool
            default = gpu ? nothing : false
        "--precond"
            help = mode == :mpi ? "Preconditioner. Defaults to cpr for blocks of mass equations or dilu otherwise." :
                "Preconditioner. Defaults to cprw where supported, cpr for other block systems, or dilu for scalar systems."
            metavar = choices(filter(x -> !(gpu && x == :jacobi) && !(mode == :mpi && x == :cprw), preconditioners))
            arg_type = Symbol
            range_tester = x -> x in preconditioners && !(gpu && x == :jacobi)
        "--global-smoother"
            help = "CPR full-system smoother. Default is dilu."
            metavar = choices(filter(!=(gpu ? :jacobi : :vendor_ilu), global_smoothers))
            arg_type = Symbol
            range_tester = x -> x in global_smoothers
        "--amg-rebuild-tolerance"
            help = "Relative pressure-matrix change needed for a scheduled full KA AMG rebuild (0 always rebuilds). Default is 0.05."
            metavar = "<fraction>"
            arg_type = Float64
            range_tester = x -> x >= 0
        "--update-type"
            help = "KA AMG reuse on full updates. Defaults to partial_sparsity with zero rebuild tolerance and memory otherwise."
            metavar = choices(amg_update_types)
            arg_type = Symbol
            range_tester = x -> x in amg_update_types
        "--update-type-partial"
            help = "KA AMG reuse on partial updates. Default is operators."
            metavar = choices(amg_update_types)
            arg_type = Symbol
            range_tester = x -> x in amg_update_types
        "--update-interval"
            help = "Schedule for full CPR updates. Default is ministep."
            metavar = choices(cpr_update_intervals)
            arg_type = Symbol
            range_tester = x -> x in cpr_update_intervals
        "--update-interval-partial"
            help = "Schedule for partial CPR updates. Default is iteration."
            metavar = choices(cpr_update_intervals)
            arg_type = Symbol
            range_tester = x -> x in cpr_update_intervals
        "--partial-update"
            help = "Allow partial CPR updates between full updates. Default is true."
            metavar = "{true,false}"
            arg_type = Bool
        "--ka-workgroupsize"
            help = "KernelAbstractions workgroup size. Default is 256."
            metavar = "<size>"
            arg_type = Int
            range_tester = x -> x > 0
        "--rtol"
            help = "Relative linear solver tolerance. Defaults to 0.005 for CPR or CPRW and 0.01 for standalone smoothers."
            metavar = "<tolerance>"
            arg_type = Float64
            range_tester = x -> 0 < x < 1
        "--tol-cnv"
            help = "CNV tolerance for mass balance. Defaults to deck TUNING or $(input_default(:tol_cnv))."
            metavar = "<tolerance>"
            arg_type = Float64
            range_tester = x -> isfinite(x) && x > 0
        "--tol-mb"
            help = "Scaled L1 mass balance tolerance. Defaults to deck TUNING or $(input_default(:tol_mb))."
            metavar = "<tolerance>"
            arg_type = Float64
            range_tester = x -> isfinite(x) && x > 0
        "--tol-dp-well"
            help = "Well connection pressure drop tolerance in bar. Default is 0.001."
            metavar = "<bar>"
            arg_type = Float64
            range_tester = x -> isfinite(x) && x > 0
        "--tol-factor-final-iteration"
            help = "Factor applied to tolerances on the final Newton iteration. Default is 1.0."
            metavar = "<factor>"
            arg_type = Float64
            range_tester = x -> isfinite(x) && x > 0
        "--relaxation"
            help = "Use nonlinear solver damping. Default is false."
            metavar = "{true,false}"
            arg_type = Bool
        "--max-nonlinear-iterations"
            help = "Maximum Newton iterations before cutting a timestep. Default is 15."
            metavar = "<iterations>"
            arg_type = Int
            range_tester = x -> 1 <= x <= 10000
        "--max-timestep-cuts"
            help = "Maximum cuts per internal timestep. Default is 25."
            metavar = "<cuts>"
            arg_type = Int
            range_tester = x -> 0 <= x <= 10000
        "--presolve-wells"
            help = "Solve wells before each linearization. Default is false."
            metavar = "{true,false}"
            arg_type = Bool
    end

    add_arg_group!(settings, "Timestepping")
    @add_arg_table! settings begin
        "--number-of-steps"
            help = "Limit the number of report steps. Runs all report steps by default."
            metavar = "<steps>"
            arg_type = Int
            range_tester = x -> x > 0
        "--timesteps"
            help = "Adaptive timestep selectors. Default is auto. Step limits and cuts still apply."
            metavar = choices(timestep_selectors)
            arg_type = Symbol
            range_tester = x -> x in timestep_selectors
        "--timesteps-target-iterations"
            help = "Target Newton iterations for timestep selection. Default is 8."
            metavar = "<iterations>"
            arg_type = Int
            range_tester = x -> x > 0
        "--timesteps-target-ds"
            help = "Target saturation change as a fraction (requires jdrun-hypre and auto timestepping). Disabled by default."
            metavar = "<fraction>"
            arg_type = Float64
            range_tester = x -> 0 < x <= 1
        "--initial-timestep"
            help = "Initial timestep in days. Defaults to deck TUNING or $(input_default(:initial_dt) / si_unit(:day)) day."
            metavar = "<days>"
            arg_type = Float64
            range_tester = x -> isfinite(x) && x > 0
        "--max-timestep"
            help = "Maximum timestep in days. Defaults to deck TUNING or $(input_default(:max_timestep) / si_unit(:day)) days."
            metavar = "<days>"
            arg_type = Float64
            range_tester = x -> x > 0
    end

    add_arg_group!(settings, "Input setup")
    if gpu
        @add_arg_table! settings begin
            "--wells-on-device"
                help = "Merge similar wells and evaluate them on the GPU"
                metavar = "{true,false}"
                arg_type = Bool
                default = true
        end
    end
    @add_arg_table! settings begin
        "--skip-tracers"
            help = "Skip tracer transport in Eclipse .DATA decks"
            action = :store_true
    end

    add_arg_group!(settings, "Output and printing")
    @add_arg_table! settings begin
        "--output-path"
            help = "Output directory (defaults to <case>_jutul_output next to the input file)"
            metavar = "<directory>"
        "--vtk"
            help = "Also export reservoir states as a ParaView VTK time series"
            action = :store_true
        "--csv"
            help = "Write field, well and group summary tables in CSV/ in SI units"
            action = :store_true
        "--output-substates"
            help = "Include internal timesteps in restart, summary, CSV and VTK output (single process)"
            action = :store_true
        "--header"
            help = "Print the JutulDarcy logo and project links at startup"
            metavar = "{true,false}"
            arg_type = Bool
            default = true
        "--info-level"
            help = "Simulation verbosity. Fractional levels are accepted. Use -1 for quiet."
            metavar = "<level>"
            arg_type = Float64
            default = 1.0
        "--report-level"
            help = "Internal Jutul report detail (0 for compact statistics, 1 or higher for iteration details)"
            metavar = "<level>"
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
    return settings
end

function parse_commandline(argv, command, mode)
    return parse_args(argv, commandline_settings(command, mode); as_symbols = true)
end

function solver_smoother(args, mode)
    precond = args[:precond]
    smoother = args[:global_smoother]
    if precond in STANDALONE_PRECONDITIONERS
        isnothing(smoother) || throw(ArgumentError(
            "--global-smoother applies to CPR. Use --precond alone to select a standalone smoother."))
        smoother = precond
    end
    gpu = mode in (:ka_cuda, :ka_amd)
    smoother == :vendor_ilu && !gpu && throw(ArgumentError("--global-smoother vendor_ilu requires a GPU launcher"))
    smoother == :jacobi && gpu && throw(ArgumentError("Jacobi is available on CPU backends"))
    smoother = something(smoother, :dilu)
    if smoother == :dilu
        return :ka_dilu
    elseif smoother == :jacobi && mode == :ka
        # KA CPU uses host vectors, but its symbolic smoother factory has no Jacobi entry.
        return JacobiPreconditioner()
    elseif smoother == :ilu
        return :ilu0
    end
    return smoother
end

function simulation_kwargs(args, mode)
    method = get(args, :method, :newton)
    if !isnothing(args[:timesteps_target_ds])
        mode == :default || throw(ArgumentError("--timesteps-target-ds requires jdrun-hypre"))
        something(args[:timesteps], :auto) == :auto ||
            throw(ArgumentError("--timesteps-target-ds requires --timesteps auto"))
    end
    if !isnothing(args[:timesteps_target_iterations]) && args[:timesteps] == :none
        throw(ArgumentError("--timesteps-target-iterations requires --timesteps auto or iteration"))
    end
    if mode == :mpi
        args[:precond] == :cprw && throw(ArgumentError("MPI supports --precond cpr for pressure preconditioning"))
        args[:output_substates] && throw(ArgumentError("MPI output supports report steps only"))
    end
    kwarg = Dict{Symbol, Any}(
        :mode => mode,
        :method => method,
        :info_level => args[:info_level],
        :report_level => args[:report_level],
        :output_substates => args[:output_substates],
        :error_on_incomplete => true,
    )
    for key in (:linear_solver, :precond, :rtol, :tol_cnv, :tol_mb, :tol_dp_well,
            :tol_factor_final_iteration, :relaxation, :max_nonlinear_iterations,
            :max_timestep_cuts, :presolve_wells, :timesteps, :ka_workgroupsize)
        value = args[key]
        isnothing(value) || (kwarg[key] = value)
    end
    if mode != :mpi
        get!(kwarg, :precond, :cprw)
    end
    solver = args[:linear_solver]
    gpu = mode in (:ka_cuda, :ka_amd)
    if gpu
        kwarg[:wells_on_device] = args[:wells_on_device]
    end
    fgmres_compatible = isnothing(solver) || solver == :fgmres
    single_precision = something(args[:single_precision_linear_solver], gpu && fgmres_compatible)
    if single_precision
        (gpu || mode == :ka) || throw(ArgumentError(
            "Single precision linear solving requires a KernelAbstractions backend"))
        fgmres_compatible || throw(ArgumentError("Single precision linear solving requires FGMRES"))
        kwarg[:single_precision_linear_solver] = true
    end
    for (key, value, scale) in (
            (:target_its, args[:timesteps_target_iterations], 1),
            (:target_ds, args[:timesteps_target_ds], 1),
            (:initial_dt, args[:initial_timestep], si_unit(:day)),
            (:max_timestep, args[:max_timestep], si_unit(:day)))
        isnothing(value) || (kwarg[key] = value * scale)
    end
    smoother_type = solver_smoother(args, mode)
    linear_solver_arg = (; smoother_type)
    if mode in (:default, :mpi)
        kwarg[:linear_solver_backend] = :cpu
        linear_solver_arg = merge(linear_solver_arg, (; amg_type = :hypre))
    end
    kwarg[:linear_solver_arg] = (; linear_solver_arg...,
        (key => args[key] for key in CPR_REBUILD_OPTIONS if !isnothing(args[key]))...)
    return kwarg
end

function run_cli(argv, command, mode; check_device = nothing, run_case = run_serial_case)
    args = parse_commandline(argv, command, mode)
    try
        path = abspath(args[:filename])
        isfile(path) || throw(ArgumentError("Input file does not exist: $path"))
        name, extension = splitext(basename(path))
        extension = lowercase(extension)
        extension in (".data", ".afi", ".jl") ||
            throw(ArgumentError("Expected a .DATA, .AFI, or .jl input file, got $path"))

        kwarg = simulation_kwargs(args, mode)
        !args[:skip_tracers] || extension == ".data" ||
            throw(ArgumentError("--skip-tracers applies to Eclipse .DATA input"))

        output_dir = args[:output_path]
        output_dir = isnothing(output_dir) ? joinpath(dirname(path), name * "_jutul_output") : abspath(output_dir)
        isnothing(check_device) || check_device()
        run_case(path, extension, output_dir, name, mode, args; kwarg)
    catch err
        if err isa ArgumentError && !args[:verbose]
            println(stderr, command, ": ", err.msg)
            return 1
        end
        rethrow()
    end
    return 0
end
