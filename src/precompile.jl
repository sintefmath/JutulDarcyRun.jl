using JutulDarcy: get_1d_reservoir, setup_well, setup_reservoir_model,
    setup_reservoir_state, setup_reservoir_forces, MultiPhaseCompositionalSystemLV,
    LiquidPhase, VaporPhase, InjectorControl, ProducerControl, TotalRateTarget,
    BottomHolePressureTarget
using JutulDarcy.MultiComponentFlash: GenericCubicEOS, MolecularProperty, MultiComponentMixture, SoreideWhitson

const PRECOMPILE_CASES = (
    (:single_phase, false),
    (:immiscible_2ph, false),
    (:deck_spe1, false),
    (:deck_spe9, false),
    (:deck_olympus, false),
    (:geothermal, true),
    (:deck_thermal, true),
    (:co2brine_immiscible, false),
    (:co2brine_immiscible, true),
    (:co2brine_kvalue, false),
    (:co2brine_kvalue, true),
    (:gas_water, false),
    (:gas_water, true),
    (:deck_tracer, false),
)

function precompile_case(physics, thermal)
    if physics in (:co2brine_immiscible, :co2brine_kvalue, :gas_water)
        return precompile_storage_case(physics, thermal)
    elseif physics in (:deck_spe1, :deck_spe9)
        name = physics == :deck_spe1 ? "SPE1" : "SPE9"
        deck = JutulDarcy.GeoEnergyIO.test_input_file_path(name, name * ".DATA")
    elseif physics == :deck_olympus
        deck = joinpath(pkgdir(@__MODULE__), "test", "data", "OLYMPUS_SMALL.DATA")
    elseif physics in (:deck_thermal, :deck_tracer)
        filename = physics == :deck_thermal ? "HEAT_SMALL.DATA" : "TRACER_SMALL.DATA"
        deck = joinpath(pkgdir(JutulDarcy), "test", "data", "thermal", filename)
    else
        return JutulDarcy.setup_mini_wellcase(physics; nstep = 1, total_time = 60.0, block_backend = true)
    end
    return load_precompile_deck(deck)[1:1]
end

function load_precompile_deck(deck; split_wells = false)
    # Keep case setup diagnostics quiet while letting failures reach the workload handler.
    return redirect_stdio(; stdout = devnull, stderr = devnull) do
        load_case(deck, ".data"; verbose = false, split_wells)
    end
end

function precompile_olympus_setup(; split_wells = false, strict = true)
    try
        deck = JutulDarcy.GeoEnergyIO.test_input_file_path("OLYMPUS_1", "OLYMPUS_1.DATA")
        case = load_precompile_deck(deck; split_wells)
        return normalize_reservoir_variables!(case)
    catch err
        err isa InterruptException && rethrow()
        @warn "OLYMPUS_1 case setup precompile workload failed" split_wells exception = (err, catch_backtrace())
        strict && rethrow()
        return nothing
    end
end

function precompile_compositional_case(ncomponents; split_wells = false)
    domain = get_1d_reservoir(4)
    wells = [setup_well(domain, [1]; name = :Injector),
        setup_well(domain, [4]; name = :Producer)]
    species = ("Methane", "n-Decane", "Ethane", "n-Propane", "n-Butane", "n-Pentane",
        "n-Hexane", "n-Heptane", "n-Octane", "n-Nonane", "n-Undecane", "n-Dodecane")
    properties = [MolecularProperty(species[mod1(i, length(species))]) for i in 1:ncomponents]
    mixture = MultiComponentMixture(properties)
    system = MultiPhaseCompositionalSystemLV(GenericCubicEOS(mixture), (LiquidPhase(), VaporPhase());
        reference_densities = (750.0, 1.0))
    model, parameters = setup_reservoir_model(domain, system;
        wells, extra_out = true, block_backend = true, split_wells)
    state0 = setup_reservoir_state(model; Pressure = 8e6, Temperature = 320.0,
        OverallMoleFractions = fill(1 / ncomponents, ncomponents))
    injection = zeros(ncomponents)
    injection[1] = 1.0
    controls = Dict(
        :Injector => InjectorControl(TotalRateTarget(1e-5), injection; density = 1.0),
        :Producer => ProducerControl(BottomHolePressureTarget(7.95e6)),
    )
    forces = setup_reservoir_forces(model; control = controls)
    return JutulCase(model, [1.0], forces; state0, parameters)
end

function precompile_storage_case(physics, thermal)
    domain = get_1d_reservoir(4)
    wells = [setup_well(domain, [1]; name = :Injector),
        setup_well(domain, [4]; name = :Producer)]
    if physics == :gas_water
        mixture = MultiComponentMixture(["Water", "Methane", "CarbonDioxide"];
            names = ["H2O", "C1", "CO2"])
        eos = GenericCubicEOS(mixture, SoreideWhitson(mixture))
        system = MultiPhaseCompositionalSystemLV(eos, (LiquidPhase(), VaporPhase());
            reference_densities = (1000.0, 1.0))
        options = (;)
        composition = (OverallMoleFractions = [0.8, 0.15, 0.05],)
        injection = [0.0, 1.0, 0.0]
    else
        immiscible = physics == :co2brine_immiscible
        system = :co2brine
        options = (co2_physics = immiscible ? :immiscible : :kvalue, co2_source = :csp11)
        composition = immiscible ? (Saturations = [0.8, 0.2],) : (OverallMoleFractions = [0.8, 0.2],)
        injection = [0.0, 1.0]
    end
    model, parameters = setup_reservoir_model(domain, system;
        wells, extra_out = true, block_backend = true, thermal, options...)
    state0 = setup_reservoir_state(model; Pressure = 1.5e7, Temperature = 325.0, composition...)
    controls = Dict(
        :Injector => InjectorControl(TotalRateTarget(1e-5), injection; density = 100.0, temperature = 320.0),
        :Producer => ProducerControl(BottomHolePressureTarget(1.49e7)),
    )
    forces = setup_reservoir_forces(model; control = controls)
    return JutulCase(model, [1.0], forces; state0, parameters)
end

function precompile_smoothers(::Val{N}, mode) where N
    # A block tridiagonal matrix exercises setup, reuse and application of the
    # full-system smoothers without compiling entire simulations for each one.
    gpu = mode in (:ka_cuda, :ka_amd)
    T = gpu ? Float32 : Float64
    diagonal = T(4) * one(Jutul.SMatrix{N, N, T})
    off_diagonal = -diagonal / T(4)
    A = Jutul.static_sparsity_sparse([1, 1, 2, 2, 2, 3, 3], [1, 2, 1, 2, 3, 2, 3],
        [diagonal, off_diagonal, off_diagonal, diagonal, off_diagonal, off_diagonal, diagonal])
    rhs = fill(one(T), 3 * N)
    if mode == :default
        context = ParallelCSRContext()
    else
        ka_backend = JutulDarcy.kernel_abstractions_backend(mode)
        context = Jutul.KernelAbstractionsContext(ka_backend; float_type = T, index_type = gpu ? Int32 : Int)
        A = Jutul.Adapt.adapt(context, A)
        rhs = Jutul.transfer(context, rhs)
    end
    backend = mode == :default ? :cpu : :ka
    array_type = mode == :ka_cuda ? "CUDA" : mode == :ka_amd ? "AMDGPU" : "Julia"
    residual = Jutul.unsafe_reinterpret(Jutul.SVector{N, T}, rhs, 3)
    output = similar(rhs)
    for name in (:ilu0, :ka_dilu, :spai0)
        smoother = JutulDarcy.reservoir_system_smoother(name; backend, array_type)
        Jutul.update_preconditioner!(smoother, A, residual, context, nothing)
        Jutul.update_preconditioner!(smoother, A, residual, context, nothing)
        Jutul.apply!(output, smoother, rhs)
    end
    return nothing
end

function precompile_simulations(mode, command; strict = true)
    precompile_olympus_setup(; strict)
    mktempdir() do dir
        args = parse_commandline(["case.jl", "--info-level", "-1"], command, mode)
        for (physics, thermal) in PRECOMPILE_CASES
            try
                case = precompile_case(physics, thermal)
                run_simulation(case, dir, "case", mode, args)
            catch err
                err isa InterruptException && rethrow()
                @warn "Simulation precompile workload failed" mode physics thermal exception = (err, catch_backtrace())
                strict && rethrow()
            end
        end
        # A one-component system uses the single-phase workload above.
        for ncomponents in 2:MAX_PRECOMPILE_COMPONENTS
            try
                case = precompile_compositional_case(ncomponents)
                run_simulation(case, dir, "case", mode, args)
            catch err
                err isa InterruptException && rethrow()
                @warn "Compositional precompile workload failed" mode ncomponents exception = (err, catch_backtrace())
                strict && rethrow()
            end
        end
        for ncomponents in 1:MAX_PRECOMPILE_COMPONENTS
            try
                precompile_smoothers(Val(ncomponents), mode)
            catch err
                err isa InterruptException && rethrow()
                @warn "Smoother precompile workload failed" mode ncomponents exception = (err, catch_backtrace())
                strict && rethrow()
            end
        end
    end
    return nothing
end
