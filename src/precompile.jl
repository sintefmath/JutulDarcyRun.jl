using JutulDarcy: get_1d_reservoir, setup_well, setup_reservoir_model,
    setup_reservoir_state, setup_reservoir_forces, MultiPhaseCompositionalSystemLV,
    LiquidPhase, VaporPhase, InjectorControl, ProducerControl, TotalRateTarget,
    BottomHolePressureTarget
using JutulDarcy.MultiComponentFlash: GenericCubicEOS, MultiComponentMixture, SoreideWhitson

function can_run_mpi_precompile()
    MPI.Finalized() && return false
    MPI.Initialized() && return MPI.Comm_size(MPI.COMM_WORLD) == 1
    # An uninitialized compiler worker can inherit a rank's launcher environment.
    return !any(key -> haskey(ENV, key), (
        "OMPI_COMM_WORLD_RANK", "PMI_RANK", "PMIX_RANK", "MPI_LOCALRANKID",
        "MV2_COMM_WORLD_RANK", "SLURM_PROCID",
    ))
end

const PRECOMPILE_CASES = (
    (:single_phase, false),
    (:immiscible_2ph, false),
    (:bo_spe1, false),
    (:compositional_2ph_3c, false),
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
    elseif physics in (:deck_thermal, :deck_tracer)
        filename = physics == :deck_thermal ? "HEAT_SMALL.DATA" : "TRACER_SMALL.DATA"
        deck = joinpath(pkgdir(JutulDarcy), "test", "data", "thermal", filename)
        return load_case(deck, ".data"; verbose = false, split_wells = false)[1:1]
    end
    return JutulDarcy.setup_mini_wellcase(physics; nstep = 1, total_time = 60.0, block_backend = true)
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

function precompile_simulations(mode; strict = true)
    mktempdir() do dir
        args = parse_commandline(["case.jl", "--info-level", "-1"], "jdrun",
            mode in (:ka, :default) ? :cpu : mode)
        for (physics, thermal) in PRECOMPILE_CASES
            args[:precond] = physics == :single_phase ? :ilu0 : nothing
            try
                case = precompile_case(physics, thermal)
                run_simulation(case, dir, "case", mode, args)
            catch err
                err isa InterruptException && rethrow()
                @warn "Simulation precompile workload failed" mode physics thermal exception = (err, catch_backtrace())
                strict && rethrow()
            end
        end
    end
    return nothing
end
