function run_mpi_case(path, extension, output_dir, name, _mode, args; kwarg)
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    function run_rank()
        case = load_case(path, extension; verbose = rank == 0 && args[:verbose],
            split_wells = true, include_tracers = !args[:skip_tracers])
        # A case script can define new methods needed for simulator setup and output.
        Base.invokelatest(run_mpi_simulation, case, output_dir, name, comm, args; kwarg)
    end
    try
        if rank == 0
            with_run_log(run_rank, output_dir, name, path)
        else
            run_rank()
        end
    catch err
        if MPI.Comm_size(comm) > 1
            failure = args[:verbose] ? CapturedException(err, catch_backtrace()) : err
            println(stderr, "jdrun-mpi rank $rank: $(sprint(showerror, failure))")
            flush(stderr)
            # A peer may already be waiting in a collective.
            MPI.Abort(comm, 1)
        end
        rethrow()
    end
    return nothing
end

function run_mpi_simulation(case::JutulCase, output_dir, name, comm, args; kwarg = simulation_kwargs(args, :mpi))
    case = limit_case(case, args[:number_of_steps])
    check_solver_layout(case, :mpi, args)
    normalize_reservoir_variables!(case)
    reservoir = JutulDarcy.reservoir_model(case)
    components = JutulDarcy.number_of_components(reservoir.system)
    block_size = degrees_of_freedom_per_entity(reservoir, Cells())
    # The distributed CPR pressure weights assume a block of mass equations.
    cpr_compatible = components > 1 && block_size == components
    precond = get(kwarg, :precond, cpr_compatible ? :cpr : :dilu)
    if !cpr_compatible && precond in (:cpr, :cprw)
        throw(ArgumentError("MPI CPR requires a block with multiple mass equations only. Use --precond dilu for this case."))
    end
    is_root = MPI.Comm_rank(comm) == 0
    scratch = if is_root
        mkpath(output_dir)
        dir = mktempdir(output_dir)
        mkdir(joinpath(dir, "jutul"))
        dir
    else
        nothing
    end
    scratch = MPI.bcast(scratch, 0, comm)
    sim, config = setup_reservoir_simulator(case; output_path = joinpath(scratch, "jutul"), kwarg..., precond)
    # Submodel construction recreates the output list from the flow system.
    local_reservoir = JutulDarcy.reservoir_model(sim.storage.model)
    union!(local_reservoir.output_variables, reservoir.output_variables, keys(reservoir.primary_variables))
    JutulDarcy.CaseValidation.validate(case; info_level = config[:info_level] - 1)
    raw = simulate!(sim, case.dt; forces = case.forces, config, start_date = case.start_date)
    # Jutul consolidates the rank-local states. All ranks may still be reading reports.
    MPI.Barrier(comm)
    if is_root
        try
            result = ReservoirSimResult(case, raw; start_date = case.start_date)
            publish_output(result, case, scratch, output_dir, name, args)
        finally
            rm(scratch; recursive = true)
        end
    end
    return nothing
end
