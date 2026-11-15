function require_backend(mode)
    if mode == :ka_cuda
        @eval import CUDA
        Base.invokelatest(() -> CUDA.functional() && CUDA.has_cuda_gpu()) ||
            throw(ArgumentError("No functional CUDA device found. jdrun-cuda requires a CUDA-capable GPU and driver."))
    elseif mode == :ka_amd
        @eval import AMDGPU
        Base.invokelatest(() -> AMDGPU.functional() && AMDGPU.has_rocm_gpu()) ||
            throw(ArgumentError("No functional ROCm device found. jdrun-amdgpu requires a ROCm-compatible GPU and driver."))
    end
    return nothing
end

function load_case(path, extension; verbose, split_wells)
    if extension == ".data"
        case = setup_case_from_data_file(path; verbose, split_wells, parse_arg = (; verbose))
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

function run_simulation(case::JutulCase, output_dir, name, mode, args)
    isempty(case.dt) && throw(ArgumentError("The case has no report steps"))
    limit = args[:number_of_steps]
    if !isnothing(limit) && limit < length(case.dt)
        case = case[1:limit]
    end
    kwarg = simulation_kwargs(args, mode)
    mkpath(output_dir)
    mktempdir(output_dir) do scratch
        native = mkdir(joinpath(scratch, "jutul"))
        result = simulate_reservoir(case; output_path = native, kwarg...)
        standard = mkdir(joinpath(scratch, "standard"))
        # Jutul reads the completed report steps back from its native files.
        # Write a complete set of standard files before replacing any previous output.
        write_standard_output(result, case, joinpath(standard, "CASE"); unified = true)
        files = [(joinpath(standard, file), joinpath(output_dir, name * splitext(file)[2]))
            for file in readdir(standard)]
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
    end
    if args[:info_level] >= 0
        jutul_message("Output", "Wrote unified Eclipse output to $output_dir")
    end
    return nothing
end
