function publish_output(result, case, scratch, output_dir, name, args)
    standard = mkdir(joinpath(scratch, "standard"))
    # Stage all formats before replacing any previous output.
    paths = write_standard_output(restart_at_output_times(result), case, joinpath(standard, "CASE"); unified = true)
    add_restart_wells(paths.restart, case)
    files = [(joinpath(standard, file), joinpath(output_dir, name * splitext(file)[2]))
        for file in readdir(standard)]
    formats = ["unified Eclipse"]
    if args[:csv]
        write_summary_csv(result, mkdir(joinpath(scratch, "CSV")))
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
        csv_path = joinpath(output_dir, "CSV")
        if args[:csv] && ispath(csv_path) && !isdir(csv_path)
            throw(ArgumentError("CSV output path is not a directory: $csv_path"))
        end
        for (source, path) in files
            Base.Filesystem.rename(source, path)
        end
        if args[:csv]
            replace_csv_directory(scratch, csv_path)
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

function restart_at_output_times(result)
    raw = result.result
    length(raw.states) == length(result.states) && return result

    # StandardIO reads restart states and times from the raw Jutul result.
    states, dt, _ = Jutul.expand_to_ministeps(raw)
    reports = [(ministeps = ((success = true, dt = step),),) for step in dt]
    expanded = Jutul.SimResult(states, reports, raw.start_timestamp)
    return JutulDarcy.ReservoirSimResult(result.wells, result.states, result.summary, result.time, expanded, result.extra)
end

function write_summary_csv(result, folder)
    values = result.summary["VALUES"]
    field = sort!(collect(pairs(values["FIELD"])); by = first)
    write_csv_table(joinpath(folder, "FIELD.csv"), result.time, field)

    wells = values["WELLS"]
    names = sort!(collect(keys(wells)))
    keywords = sort!(collect(Set(keyword for well in names for keyword in keys(wells[well]))))
    blanks = fill("", length(result.time))
    for keyword in keywords
        columns = [well => get(wells[well], keyword, blanks) for well in names]
        write_csv_table(joinpath(folder, keyword * ".csv"), result.time, columns)
    end

    groups = ["$keyword:$group" => samples
        for (group, series) in values["GROUP"] for (keyword, samples) in series]
    if !isempty(groups)
        sort!(groups; by = first)
        write_csv_table(joinpath(folder, "GROUP.csv"), result.time, groups)
    end
    return nothing
end

function write_csv_table(path, time, columns)
    for (name, samples) in columns
        length(samples) == length(time) || throw(ArgumentError("CSV series $name in $path has a different length than the time axis"))
    end
    return open(path, "w") do io
        writedlm(io, (["time_seconds"; first.(columns)],), ',')
        # Iterate by row without materializing a matrix or a zip of every well column.
        rows = (Iterators.flatten(((t,), (last(column)[i] for column in columns))) for (i, t) in enumerate(time))
        writedlm(io, rows, ',')
    end
end

function replace_csv_directory(scratch, destination)
    previous = joinpath(scratch, "previous-csv")
    if isdir(destination)
        Base.Filesystem.rename(destination, previous)
    end
    try
        Base.Filesystem.rename(joinpath(scratch, "CSV"), destination)
    catch
        isdir(previous) && Base.Filesystem.rename(previous, destination)
        rethrow()
    end
    return nothing
end

function write_vtk_output(result, case, folder, name)
    @eval import WriteVTK
    mesh = physical_representation(JutulDarcy.reservoir_domain(case))
    files = Base.invokelatest(export_mesh_vtu, mesh, name * ".pvd"; folder, cell_data = result.states)
    collection = joinpath(folder, name * ".pvd")

    # Jutul numbers the collection entries by state index. Use elapsed seconds.
    with_elapsed_time(index) = string(result.time[Int(parse(Float64, index))])
    contents = replace(read(collection, String), r"(?<=\btimestep=\")[^\"]+" => with_elapsed_time)
    write(collection, contents)
    return files
end
