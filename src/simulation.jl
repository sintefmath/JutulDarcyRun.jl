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
    if any(key -> !isnothing(args[key]), (:linear_solver, :precond, :global_smoother, :rtol))
        throw(ArgumentError("Krylov solver options require a cell-major case. Set block_backend = true when building the model."))
    end
    return nothing
end

function run_simulation(case::JutulCase, output_dir, name, mode, args; kwarg = simulation_kwargs(args, mode))
    case = limit_case(case, args[:number_of_steps])
    check_solver_layout(case, mode, args)
    mkpath(output_dir)
    mktempdir(output_dir) do scratch
        native = mkdir(joinpath(scratch, "jutul"))
        result = simulate_reservoir(case; output_path = native, kwarg...)
        publish_output(result, case, scratch, output_dir, name, args)
    end
    return nothing
end

function publish_output(result, case, scratch, output_dir, name, args)
    standard = mkdir(joinpath(scratch, "standard"))
    # Write a complete set of standard files before replacing any previous output.
    write_standard_output(result, case, joinpath(standard, "CASE"); unified = true)
    files = [(joinpath(standard, file), joinpath(output_dir, name * splitext(file)[2]))
        for file in readdir(standard)]
    formats = ["unified Eclipse"]
    if args[:csv]
        csv = name * ".csv"
        source = joinpath(standard, csv)
        write_summary_csv(result, source)
        push!(files, (source, joinpath(output_dir, csv)))
        push!(formats, "CSV")
    end
    if args[:vtk]
        vtk_files = write_vtk_output(result, case, standard, name)
        # Publish the collection only after the files it references are in place.
        sort!(vtk_files; by = path -> endswith(path, ".pvd"))
        for source in vtk_files
            push!(files, (source, joinpath(output_dir, basename(source))))
        end
        push!(formats, "VTK")
    end
    mkpidlock(joinpath(output_dir, ".jdrun.lock"); stale_age = 60) do
        for (_, path) in files
            isdir(path) && throw(ArgumentError("Output file path is a directory: $path"))
        end
        for (source, path) in files
            Base.Filesystem.rename(source, path)
        end
    end
    if args[:print_wells]
        result.wells()
    end
    if args[:info_level] >= 0
        jutul_message("Output", "Wrote $(join(formats, ", ")) output to $output_dir")
    end
    return nothing
end

function write_summary_csv(result, path)
    values = result.summary["VALUES"]
    scopes = (
        ("field", Dict("FIELD" => values["FIELD"])),
        ("group", values["GROUP"]),
        ("well", values["WELLS"]),
    )
    rows = ((time, scope, name, keyword, value)
        for (scope, entities) in scopes
        for (name, series) in sort!(collect(pairs(entities)); by = first)
        for (keyword, samples) in sort!(collect(pairs(series)); by = first)
        for (time, value) in zip(result.time, samples))
    return open(path, "w") do io
        println(io, "time_seconds,scope,name,keyword,value")
        writedlm(io, rows, ',')
    end
end

function write_vtk_output(result, case, folder, name)
    @eval import WriteVTK
    mesh = physical_representation(JutulDarcy.reservoir_domain(case))
    return Base.invokelatest(export_mesh_vtu, mesh, name * ".pvd";
        folder, cell_data = result.states)
end
