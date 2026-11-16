const STANDALONE_PRECONDITIONERS = (:ilu0, :dilu, :jacobi, :spai0)
const GLOBAL_SMOOTHERS = (:ilu0, :ilu, :dilu, :jacobi, :spai0, :vendor_ilu, :ka_ilu0, :ka_dilu, :ka_spai0)
const AMG_UPDATE_TYPES = (:memory, :sparsity, :partial_sparsity, :partial_operators, :operators)
const CPR_UPDATE_INTERVALS = (:once, :iteration, :ministep, :step)
const CPR_REBUILD_OPTIONS = (:amg_rebuild_tolerance, :update_type, :update_type_partial,
    :update_interval, :update_interval_partial, :partial_update)
const DEFAULT_GLOBAL_SMOOTHER = :dilu

function commandline_settings(command, mode)
    gpu = mode in (:ka_cuda, :ka_amd)
    precond_help = if gpu
        "Preconditioner (cprw, cpr, ilu0, dilu or spai0). Defaults to cprw where supported, cpr for other block systems, or dilu for scalar systems."
    elseif mode == :mpi
        "Preconditioner (cpr, ilu0, dilu, jacobi or spai0). Defaults to cpr for blocks of mass equations or dilu otherwise."
    else
        "Preconditioner (cprw, cpr, ilu0, dilu, jacobi or spai0). Defaults to cprw where supported, cpr for other block systems, or dilu for scalar systems."
    end
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
                help = "Nonlinear method (newton or nldd)"
                arg_type = Symbol
                default = :newton
                range_tester = x -> x in (:newton, :nldd)
        end
    end
    @add_arg_table! settings begin
        "--linear-solver"
            help = gpu ? "Krylov solver (fgmres by default, gmres or bicgstab)" :
                "Krylov solver (gmres by default, fgmres or bicgstab)"
            arg_type = Symbol
            range_tester = x -> x in (:gmres, :fgmres, :bicgstab)
        "--single-precision-linear-solver"
            help = gpu ? "Use Float32 (true or false). Defaults to true with FGMRES and false with other solvers." :
                "Use Float32 with FGMRES on CPU KA (true or false). Default is false."
            arg_type = Bool
        "--precond"
            help = precond_help
            arg_type = Symbol
            range_tester = x -> x in (:cprw, :cpr, STANDALONE_PRECONDITIONERS...) && !(gpu && x == :jacobi)
        "--global-smoother"
            help = gpu ? "CPR full-system smoother (dilu, ilu0, spai0 or vendor_ilu). Default is $DEFAULT_GLOBAL_SMOOTHER." :
                "CPR full-system smoother (dilu, ilu0, jacobi or spai0). Default is $DEFAULT_GLOBAL_SMOOTHER."
            arg_type = Symbol
            range_tester = x -> x in GLOBAL_SMOOTHERS
        "--amg-rebuild-tolerance"
            help = "Relative pressure-matrix change needed for a scheduled full KA AMG rebuild (0 always rebuilds). Default is 0.05."
            arg_type = Float64
            range_tester = x -> x >= 0
        "--update-type"
            help = "KA AMG reuse on full updates (memory, sparsity, partial_sparsity, partial_operators or operators). Defaults to partial_sparsity with zero rebuild tolerance and memory otherwise."
            arg_type = Symbol
            range_tester = x -> x in AMG_UPDATE_TYPES
        "--update-type-partial"
            help = "KA AMG reuse on partial updates (same choices as --update-type). Default is operators."
            arg_type = Symbol
            range_tester = x -> x in AMG_UPDATE_TYPES
        "--update-interval"
            help = "Schedule for full CPR updates (once, iteration, ministep or step). Default is ministep."
            arg_type = Symbol
            range_tester = x -> x in CPR_UPDATE_INTERVALS
        "--update-interval-partial"
            help = "Schedule for partial CPR updates (same choices as --update-interval). Default is iteration."
            arg_type = Symbol
            range_tester = x -> x in CPR_UPDATE_INTERVALS
        "--partial-update"
            help = "Allow partial CPR updates between full updates (true or false). Default is true."
            arg_type = Bool
        "--ka-workgroupsize"
            help = "KernelAbstractions workgroup size. Default is 256."
            arg_type = Int
            range_tester = x -> x > 0
        "--rtol"
            help = "Relative linear solver tolerance. Defaults to 0.005 for CPR or CPRW and 0.01 for standalone smoothers."
            arg_type = Float64
            range_tester = x -> 0 < x < 1
        "--tol-cnv"
            help = "CNV tolerance for mass balance. Defaults to deck TUNING or $(input_default(:tol_cnv))."
            arg_type = Float64
            range_tester = x -> isfinite(x) && x > 0
        "--tol-mb"
            help = "Scaled L1 mass balance tolerance. Defaults to deck TUNING or $(input_default(:tol_mb))."
            arg_type = Float64
            range_tester = x -> isfinite(x) && x > 0
        "--tol-dp-well"
            help = "Well connection pressure drop tolerance in bar. Default is 0.001."
            arg_type = Float64
            range_tester = x -> isfinite(x) && x > 0
        "--tol-factor-final-iteration"
            help = "Factor applied to tolerances on the final Newton iteration. Default is 1.0."
            arg_type = Float64
            range_tester = x -> isfinite(x) && x > 0
        "--relaxation"
            help = "Use nonlinear solver damping (true or false). Default is false."
            arg_type = Bool
        "--max-nonlinear-iterations"
            help = "Maximum Newton iterations before cutting a timestep. Default is 15."
            arg_type = Int
            range_tester = x -> 1 <= x <= 10000
        "--max-timestep-cuts"
            help = "Maximum cuts per internal timestep. Default is 25."
            arg_type = Int
            range_tester = x -> 0 <= x <= 10000
        "--presolve-wells"
            help = "Solve wells before each linearization (true or false). Default is false."
            arg_type = Bool
    end

    add_arg_group!(settings, "Timestepping")
    @add_arg_table! settings begin
        "--number-of-steps"
            help = "Limit the number of report steps. Runs all report steps by default."
            arg_type = Int
            range_tester = x -> x > 0
        "--timesteps"
            help = "Adaptive timestep selectors (auto, iteration or none). Default is auto. Step limits and cuts still apply."
            arg_type = Symbol
            range_tester = x -> x in (:auto, :iteration, :none)
        "--timesteps-target-iterations"
            help = "Target Newton iterations for timestep selection. Default is 8."
            arg_type = Int
            range_tester = x -> x > 0
        "--timesteps-target-ds"
            help = "Target saturation change as a fraction (requires jdrun-hypre and auto timestepping). Disabled by default."
            arg_type = Float64
            range_tester = x -> 0 < x <= 1
        "--initial-timestep"
            help = "Initial timestep in days. Defaults to deck TUNING or $(input_default(:initial_dt) / si_unit(:day)) day."
            arg_type = Float64
            range_tester = x -> isfinite(x) && x > 0
        "--max-timestep"
            help = "Maximum timestep in days. Defaults to deck TUNING or $(input_default(:max_timestep) / si_unit(:day)) days."
            arg_type = Float64
            range_tester = x -> x > 0
    end

    add_arg_group!(settings, "Input setup")
    if gpu
        @add_arg_table! settings begin
            "--wells-on-device"
                help = "Merge similar wells and evaluate them on the GPU"
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
        "--vtk"
            help = "Also export reservoir states as a ParaView VTK time series"
            action = :store_true
        "--csv"
            help = "Also export field, well and group summary time series as CSV in SI units"
            action = :store_true
        "--output-substates"
            help = "Include internal timesteps in summary, CSV, and VTK output (single process)"
            action = :store_true
        "--info-level"
            help = "Simulation verbosity. Fractional levels are accepted. Use -1 for quiet."
            arg_type = Float64
            default = 1.0
        "--report-level"
            help = "Internal Jutul report detail (0 for compact statistics, 1 or higher for iteration details)"
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
    smoother = something(smoother, DEFAULT_GLOBAL_SMOOTHER)
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
    if gpu && isnothing(solver)
        kwarg[:linear_solver] = :fgmres
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
