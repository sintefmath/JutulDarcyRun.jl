using ArgParse
using PackageCompiler: create_app
using Pkg
using TOML

function build_settings(variant, output)
    settings = ArgParseSettings(
        prog = "julia --project=executables executables/$variant/build.jl",
        description = "Build the executable and warm up on the first report step of each supplied case.",
        autofix_names = true,
    )
    @add_arg_table! settings begin
        "--output-path"
            help = "Destination for the executable bundle"
            default = output
        "--cpu-target"
            help = "Julia CPU target (use generic for wider CPU compatibility)"
            default = "native"
        "--include-lazy-artifacts"
            help = "Bundle lazy dependency artifacts as well, including optional datasets"
            action = :store_true
        "cases"
            help = "Paths to .DATA, .AFI, or .jl cases for precompilation"
            nargs = '*'
    end
    return settings
end

function stage_app(stage, app, project, repo)
    package = joinpath(stage, "app")
    mkdir(package)
    cp(joinpath(app, "src"), joinpath(package, "src"))
    sources = project["sources"]
    sources["JutulDarcyRun"] = Dict("path" => repo)
    for source in values(sources)
        if haskey(source, "path")
            source["path"] = abspath(app, source["path"])
        end
    end

    # Julia case files can import the modeling packages in the standalone bundle.
    runner = TOML.parsefile(joinpath(repo, "Project.toml"))
    for name in ("Jutul", "JutulDarcy"), table in ("deps", "compat")
        get!(project[table], name, runner[table][name])
    end
    open(io -> TOML.print(io, project), joinpath(package, "Project.toml"), "w")

    preferences_path = joinpath(app, "LocalPreferences.toml")
    preferences = isfile(preferences_path) ? TOML.parsefile(preferences_path) : Dict{String, Any}()
    # PackageCompiler traces the requested cases. Skip the launcher's separate Pkg workload.
    get!(Dict{String, Any}, preferences, project["name"])["precompile_workload"] = false
    open(io -> TOML.print(io, preferences), joinpath(package, "LocalPreferences.toml"), "w")
    return package
end

function write_warmup(stage, number, package_name, case_path)
    # PackageCompiler runs each file in a fresh process, including when several
    # Julia case scripts define a top-level case or methods with the same name.
    path = joinpath(stage, "warmup-$number.jl")
    open(path, "w") do io
        println(io, "import ", package_name)
        println(io, "include(", repr(joinpath(@__DIR__, "warmup.jl")), ")")
        println(io, "warmup_case(", package_name, ".main, ", repr(case_path), ")")
    end
    return path
end

function build_executable(directory, argv = ARGS)
    repo = dirname(@__DIR__)
    variant = basename(normpath(directory))
    app = joinpath(repo, "apps", variant)
    project = TOML.parsefile(joinpath(app, "Project.toml"))
    command = only(keys(project["apps"]))
    settings = build_settings(variant, joinpath(directory, "build", command))
    args = parse_args(argv, settings; as_symbols = true)
    cases = abspath.(args[:cases])
    for path in cases
        isfile(path) || throw(ArgumentError("Case file does not exist: $path"))
        lowercase(splitext(path)[2]) in (".data", ".afi", ".jl") ||
            throw(ArgumentError("Expected a .DATA, .AFI, or .jl case, got $path"))
    end
    destination = abspath(args[:output_path])
    isempty(basename(destination)) && (destination = dirname(destination))
    ispath(destination) && throw(ArgumentError("Output path already exists: $destination"))
    parent = mkpath(dirname(destination))

    withenv("JULIA_NUM_THREADS" => get(ENV, "JULIA_NUM_THREADS", "8")) do
        mktempdir() do stage
            package = stage_app(stage, app, project, repo)
            Pkg.activate(package) do
                Pkg.instantiate(; allow_autoprecomp = false)
            end
            workload = String[write_warmup(stage, i, project["name"], path) for (i, path) in enumerate(cases)]
            mktempdir(parent) do build_dir
                bundle = joinpath(build_dir, "bundle")
                create_app(package, bundle;
                    executables = [command => "julia_main"],
                    precompile_execution_file = workload,
                    incremental = true,
                    cpu_target = args[:cpu_target],
                    include_lazy_artifacts = args[:include_lazy_artifacts],
                    sysimage_build_args = `-O3`)
                mv(bundle, destination)
            end
        end
    end
    println("Built ", joinpath(destination, "bin", command * (Sys.iswindows() ? ".exe" : "")))
    return nothing
end
