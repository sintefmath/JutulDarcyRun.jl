function load_case(path, extension; verbose, split_wells, include_tracers = true)
    if extension == ".data"
        case = setup_case_from_data_file(path; verbose, split_wells, include_tracers, parse_arg = (; verbose))
    elseif extension == ".afi"
        case = setup_case_from_afi(path; verbose, split_wells)
    else
        Base.include(Main, path)
        case = Base.invokelatest() do
            isdefined(Main, :case) || throw(ArgumentError("$path must define a top-level variable named case"))
            return Main.case
        end
    end
    case isa JutulCase || throw(ArgumentError("$path must define case as a JutulCase, got $(typeof(case))"))
    return case
end

function limit_case(case::JutulCase, limit)
    isempty(case.dt) && throw(ArgumentError("The case has no report steps"))
    if !isnothing(limit) && limit < length(case.dt)
        case = case[1:limit]
    end
    return case
end

function check_solver_layout(case, mode, args)
    reservoir = JutulDarcy.reservoir_model(case)
    is_cell_major(matrix_layout(reservoir.context)) && return nothing
    if mode != :default
        throw(ArgumentError("This backend requires a cell-major case. Set block_backend = true when building the model."))
    end
    if any(key -> !isnothing(args[key]), (:linear_solver, :precond, :global_smoother, :rtol, CPR_REBUILD_OPTIONS...))
        throw(ArgumentError("Krylov solver options require a cell-major case. Set block_backend = true when building the model."))
    end
    return nothing
end

function case_preconditioner(case, args, precond)
    if isnothing(args[:precond]) && precond == :cprw
        reservoir = JutulDarcy.reservoir_model(case)
        if degrees_of_freedom_per_entity(reservoir, Cells()) == 1
            return :dilu
        end
        model = case.model
        grouped = model isa MultiModel && model.groups isa AbstractVector && maximum(model.groups; init = 0) > 1
        # The well pressure map uses off-diagonal CSR blocks and assumes no tracer equations.
        if !grouped || !(reservoir.context isa ParallelCSRContext) || haskey(reservoir.equations, :tracers)
            return :cpr
        end
    end
    return precond
end

function normalize_reservoir_variables!(case)
    reservoir = JutulDarcy.reservoir_model(case)
    multi = case.model isa MultiModel
    state = multi ? case.state0[:Reservoir] : case.state0
    if haskey(reservoir.parameters, :FluidVolume)
        parameters = multi ? case.parameters[:Reservoir] : case.parameters
        volume = parameters[:FluidVolume]
        JutulDarcy.set_rock_compressibility!(reservoir; compressibility = 0.0)
        delete!(parameters, :FluidVolume)
        parameters[:StaticFluidVolume] = volume
        state[:FluidVolume] = copy(volume)
    end
    nph = JutulDarcy.number_of_phases(reservoir.system)
    if nph > 1 && isnothing(get_variable(reservoir, :CapillaryPressure; throw = false))
        zero_pressure = get_1d_interpolator([0.0, 1.0], [0.0, 0.0])
        pc = ntuple(_ -> zero_pressure, nph - 1)
        Jutul.set_secondary_variables!(reservoir; CapillaryPressure = JutulDarcy.SimpleCapillaryPressure(pc))
        pressure = state[:Pressure]
        state[:CapillaryPressure] = zeros(eltype(pressure), nph - 1, length(pressure))
    end
    return case
end

function run_serial_case(path, extension, output_dir, name, mode, args; kwarg)
    return with_run_log(output_dir, name, path; header = args[:header]) do
        case = load_case(path, extension; verbose = args[:verbose], split_wells = kwarg[:method] == :nldd,
            include_tracers = !args[:skip_tracers])
        # Julia case scripts can define methods needed by simulator setup and output.
        Base.invokelatest(run_simulation, case, output_dir, name, mode, args; kwarg)
    end
end

function run_simulation(case::JutulCase, output_dir, name, mode, args; kwarg = simulation_kwargs(args, mode))
    case = limit_case(case, args[:number_of_steps])
    check_solver_layout(case, mode, args)
    normalize_reservoir_variables!(case)
    if get(kwarg, :wells_on_device, false) && case.model isa MultiModel
        try
            case = JutulDarcy.merge_similar_wells(case)
        catch err
            err isa ArgumentError || rethrow()
            kwarg = copy(kwarg)
            kwarg[:wells_on_device] = false
            if args[:info_level] >= 0
                jutul_message("Wells", "Evaluating wells on the CPU. $(err.msg)")
            end
        end
    end
    precond = case_preconditioner(case, args, kwarg[:precond])
    mkpath(output_dir)
    mktempdir(output_dir) do scratch
        native = mkdir(joinpath(scratch, "jutul"))
        result = simulate_reservoir(case; output_path = native, kwarg..., precond)
        publish_output(result, case, scratch, output_dir, name, args)
    end
    return nothing
end
